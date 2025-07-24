
library(data.table)   
library(rugarch)      
library(readxl)          


dt_long <- read_excel('C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/data/data.xlsx')          
dt_long <- as.data.table(dt_long)
dt_long[, Date := as.Date(Date, format = "%d/%m/%Y")]

dt_wide <- dcast(dt_long, Date ~ Ticker, value.var = "Ret")
setorder(dt_wide, Date)                           

ret_xts <- xts(dt_wide[, -1, with = FALSE], order.by = dt_wide$Date)


spec_norm <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE),  # AR(1)
  distribution.model = "norm"                                       # Normal
)

## FIT TICKER-BY-TICKER, COLLECT RESIDUALS & PITs

ticker_names  <- colnames(ret_xts)
resid_list    <- list()   # standardised residuals
pit_list      <- list()   # PITs

for (sym in ticker_names) {
  cat("Fitting", sym, "...\n")
  y <- na.omit(ret_xts[, sym])                      # drop leading/trailing NAs
  
  fit <- ugarchfit(spec = spec_norm, data = y)
  
  # standardised (filtered) residuals: ε_t / σ_t
  z_t <- residuals(fit, standardize = TRUE)
  
  # PIT under Normal assumption
  u_t <- pnorm(z_t)
  
  resid_list[[sym]] <- z_t
  pit_list[[sym]]   <- u_t
}

# Align by date and bind back into xts/data.table
resid_xts <- do.call(merge, resid_list)   # one column per ticker
pit_xts   <- do.call(merge, pit_list)

# If you prefer data.tables:
resid_dt <- as.data.table(resid_xts, keep.rownames = "Date")
pit_dt   <- as.data.table(pit_xts,   keep.rownames = "Date")

# Quick peek
head(resid_dt)
head(pit_dt)























