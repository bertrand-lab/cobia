# model fitting for estimating ion peak width

# data retrieved from Hsieh et al 2013 'Effects of Column and Gradient Lengths on Peak Capacity and Peptide Identification in Nanoflow LC-MS/MS of Complex Proteomic Samples' Supplementary Table 2

# Data then used to estimate ion peak width from LC column characteristics

peak_width <- c(0.197, 0.164, 0.138, 0.11, 0.308, 0.261, 0.209, 0.161, 0.412, 0.343, 0.272, 0.208)
column_len <- c(rep(100, 3), rep(200, 3), rep(400, 3), rep(600, 3))
gradient <- c(rep(30, 4), rep(60, 4), rep(90, 4))
df <- data.frame(column_len, gradient, peak_width)
plot(df)

lm_out1 <- lm(peak_width ~ column_len + gradient)

#kleiner 260
predict(object = lm_out1, newdata = data.frame(column_len = 500, gradient = 260))

#kleiner 460
predict(object = lm_out1, newdata = data.frame(column_len = 500, gradient = 460))

# broberg
predict(object = lm_out1, newdata = data.frame(column_len = 250, gradient = 129))

# aylward
predict(object = lm_out1, newdata = data.frame(column_len = 600, gradient = 100))

## own data

predict(object = lm_out1, newdata = data.frame(column_len = , gradient = 120))
