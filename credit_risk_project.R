library(tidyverse)
library(MASS)

# Resolve masking conflicts
select <- dplyr::select
filter <- dplyr::filter

# =========================
# 1. LOAD DATA
# =========================
data <- read.csv("data/UCI_Credit_Card.csv")
data$default <- data$default.payment.next.month

# =========================
# 2. SANITIZATION
# =========================
data[data == Inf | data == -Inf] <- NA

for (i in names(data)) {
  if (is.numeric(data[[i]])) {
    data[[i]][is.na(data[[i]])] <- median(data[[i]], na.rm = TRUE)
  }
}

data <- data[, sapply(data, function(x) !is.numeric(x) || sd(x, na.rm = TRUE) > 0)]

# =========================
# 3. SAFE SCALING
# =========================
safe_scale <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# =========================
# 4. LATENT SCORE + RATINGS
# =========================
set.seed(123)
data$score <-
  safe_scale(data$LIMIT_BAL) * 0.5 +
  safe_scale(data$AGE) * (-0.2) +
  rnorm(nrow(data), 0, 1)

data$rating <- cut(
  data$score,
  breaks = 8,
  labels = c("D","CCC","B","BB","BBB","A","AA","AAA"),
  include.lowest = TRUE
)
data$rating <- factor(data$rating, ordered = TRUE)

stopifnot(!any(is.na(data$rating)))

# =========================
# 5. ORDINAL PROBIT
# =========================
model_data <- data %>%
  dplyr::select(rating, LIMIT_BAL, AGE) %>%
  mutate(
    LIMIT_BAL = safe_scale(LIMIT_BAL),
    AGE       = safe_scale(AGE)
  ) %>%
  filter(complete.cases(.))

fit <- polr(
  rating ~ LIMIT_BAL + AGE,
  data   = model_data,
  method = "probit",
  Hess   = TRUE,
  control = list(maxit = 300)
)
print(summary(fit))

# =========================
# 6. MERTON STRUCTURAL MODEL
# =========================
# За всеки клиент симулираме стойност на активите
# V_A = Asset value, D = Debt (= LIMIT_BAL като proxy)
# PD_merton = P(V_A < D) = Phi(-d2)
# d2 = (ln(V/D) + (mu - 0.5*sigma^2)*T) / (sigma*sqrt(T))

set.seed(42)
T       <- 1          # хоризонт 1 година
mu      <- 0.05       # очаквана възвръщаемост на активите
sigma   <- 0.20       # волатилност на активите

# Proxy: V = LIMIT_BAL * (1 + uniform noise), D = 0.8 * LIMIT_BAL
data$V  <- data$LIMIT_BAL * runif(nrow(data), 1.0, 1.5)
data$D  <- data$LIMIT_BAL * 0.8

data$d2 <- (log(data$V / data$D) + (mu - 0.5 * sigma^2) * T) /
           (sigma * sqrt(T))

data$PD_merton <- pnorm(-data$d2)

cat("\n--- Merton Model ---\n")
cat(sprintf("Mean PD (Merton): %.4f\n", mean(data$PD_merton)))
cat(sprintf("Mean PD (Empirical): %.4f\n", mean(data$default)))

# =========================
# 7. EL / UL
# =========================
LGD <- 0.45
EAD <- mean(data$LIMIT_BAL)
PD  <- mean(data$default)

EL  <- PD * LGD * EAD
UL  <- sqrt(PD * (1 - PD)) * LGD * EAD

cat(sprintf("\nEL: %.2f TWD | UL: %.2f TWD\n", EL, UL))

# =========================
# 8. ASSET CORRELATION (Basel II proxy)
# =========================
# Basel II formula за retail експозиции:
# rho = 0.03 * (1 - exp(-35*PD))/(1 - exp(-35)) +
#       0.16 * (1 - (1 - exp(-35*PD))/(1 - exp(-35)))

rho <- 0.03 * (1 - exp(-35 * PD)) / (1 - exp(-35)) +
       0.16 * (1 - (1 - exp(-35 * PD)) / (1 - exp(-35)))

cat(sprintf("Asset Correlation (Basel II retail): %.4f\n", rho))

# =========================
# 9. MONTE CARLO С ASSET CORRELATION
# =========================
set.seed(123)
n_loans <- 1000
n_sim   <- 5000

# Vasicek single-factor model:
# Default if: sqrt(rho)*Z + sqrt(1-rho)*epsilon < Phi^{-1}(PD)
threshold <- qnorm(PD)

losses <- replicate(n_sim, {
  Z       <- rnorm(1)                          # systematic factor
  epsilon <- rnorm(n_loans)                    # idiosyncratic
  asset_returns <- sqrt(rho) * Z + sqrt(1 - rho) * epsilon
  defaults <- as.integer(asset_returns < threshold)
  sum(defaults * LGD * EAD)
})

loss_df <- data.frame(loss = losses)

VaR_95  <- quantile(losses, 0.95)
VaR_99  <- quantile(losses, 0.99)
CVaR_99 <- mean(losses[losses >= VaR_99])

cat(sprintf("VaR 95%%: %.2f | VaR 99%%: %.2f | CVaR 99%%: %.2f\n",
            VaR_95, VaR_99, CVaR_99))

# =========================
# GRAPH 1 — PD по рейтинги
# =========================
pd_tbl <- data %>%
  group_by(rating) %>%
  summarise(PD = mean(default), N = n(), .groups = "drop")

g1 <- ggplot(pd_tbl, aes(x = rating, y = PD, fill = rating)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = scales::percent(PD, accuracy = 0.1)),
            vjust = -0.5, size = 3.5) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = "Вероятност за дефолт по кредитен рейтинг",
    subtitle = "Ordinal Probit модел — UCI Credit Card Dataset",
    x = "Рейтинг", y = "PD (%)"
  ) +
  theme_minimal(base_size = 13)

ggsave("graph1_pd_by_rating.png", g1, width = 9, height = 5)

# =========================
# GRAPH 2 — EL vs UL
# =========================
risk <- data.frame(
  Metric = c("Expected Loss (EL)", "Unexpected Loss (UL)"),
  Value  = c(EL, UL)
)

g2 <- ggplot(risk, aes(x = Metric, y = Value, fill = Metric)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(round(Value, 0))),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Expected Loss vs Unexpected Loss",
    subtitle = sprintf("LGD = %.0f%% | EAD = %s TWD | PD = %.2f%%",
                       LGD * 100, scales::comma(round(EAD)), PD * 100),
    x = "", y = "Стойност (TWD)"
  ) +
  theme_minimal(base_size = 13)

ggsave("graph2_el_ul.png", g2, width = 7, height = 5)

# =========================
# GRAPH 3 — Loss Distribution (Vasicek)
# =========================
g3 <- ggplot(loss_df, aes(x = loss)) +
  geom_histogram(bins = 50, fill = "gray60", color = "white", alpha = 0.85) +
  geom_vline(xintercept = EL,      color = "green4", linetype = "solid",  linewidth = 1.1) +
  geom_vline(xintercept = VaR_95,  color = "orange", linetype = "dashed", linewidth = 1.1) +
  geom_vline(xintercept = VaR_99,  color = "red",    linetype = "dashed", linewidth = 1.1) +
  annotate("text", x = EL,     y = Inf, label = "EL",      vjust = 2, hjust = -0.2, color = "green4") +
  annotate("text", x = VaR_95, y = Inf, label = "VaR 95%", vjust = 2, hjust = -0.2, color = "orange") +
  annotate("text", x = VaR_99, y = Inf, label = "VaR 99%", vjust = 2, hjust = -0.2, color = "red") +
  labs(
    title    = "Monte Carlo загуби — Vasicek single-factor модел",
    subtitle = sprintf("ρ (asset correlation) = %.4f | n = %d симулации", rho, n_sim),
    x = "Портфейлна загуба (TWD)", y = "Честота"
  ) +
  theme_minimal(base_size = 13)

ggsave("graph3_loss_distribution.png", g3, width = 9, height = 5)

# =========================
# GRAPH 4 — Merton PD разпределение
# =========================
g4 <- ggplot(data, aes(x = PD_merton)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.85) +
  geom_vline(xintercept = mean(data$PD_merton),
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = mean(data$PD_merton), y = Inf,
           label = sprintf("Средна PD = %.3f", mean(data$PD_merton)),
           vjust = 2, hjust = -0.1, color = "red") +
  labs(
    title    = "Разпределение на PD — Merton структурен модел",
    subtitle = sprintf("σ = %.2f | μ = %.2f | T = %d год.", sigma, mu, T),
    x = "Вероятност за дефолт (Merton)", y = "Брой клиенти"
  ) +
  theme_minimal(base_size = 13)

ggsave("graph4_merton_pd.png", g4, width = 9, height = 5)


# =========================
# ЗАПИС НА РЕЗУЛТАТИ В results/
# =========================
dir.create("results", showWarnings = FALSE)

# PD по рейтинги
write.csv(pd_tbl, "results/pd_by_rating.csv", row.names = FALSE)

# EL / UL стойности
write.csv(risk, "results/el_ul_values.csv", row.names = FALSE)

# Loss distribution статистики
loss_stats <- data.frame(
  Metric = c("EL", "VaR_95", "VaR_99", "CVaR_99", "Asset_Correlation"),
  Value  = c(EL, VaR_95, VaR_99, CVaR_99, rho)
)
write.csv(loss_stats, "results/loss_distribution_stats.csv", row.names = FALSE)

# Графики в results/
ggsave("results/graph1_pd_by_rating.png",    g1, width = 9, height = 5)
ggsave("results/graph2_el_ul.png",           g2, width = 7, height = 5)
ggsave("results/graph3_loss_distribution.png", g3, width = 9, height = 5)
ggsave("results/graph4_merton_pd.png",       g4, width = 9, height = 5)

# SessionInfo
sink("results/SessionInfo.txt")
sessionInfo()
sink()

cat("\n✅ Всички резултати са записани в папка results/\n")