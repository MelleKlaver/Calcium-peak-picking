
library(ggplot2)
library(dplyr)

#load data
data = read.csv("sCaE_example.csv")
data = data[c(0:35001),]


#y = data$Trace..1..mV.
#x = data$Time.from.zero..s.
y=data[,2]
x=data[,1]
cellname="2230MT_IL6_Ca"
original_x = data[,1]
original_y = data[,2]
  
peakdet = function(x, y, w=1, ...) {
  require(zoo)
  n = length(y)
  y.smooth = loess(y ~ x, ...)$fitted
  y.max = rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta = y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max = which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

peakdet_test = function(w, span, filename) {
  peaks = peakdet(x, y, w=w, span=span)
  
  plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
  lines(x, peaks$y.hat,  lwd=2)
  y.min = min(y)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]), col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
  
  # Save plot
  png(filename, width = 8, height = 5, units = "in", res = 600)
  plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
  lines(x, peaks$y.hat,  lwd=2)
  y.min = min(y)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]), col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
  dev.off()
}

p = ggplot(data = data.frame(x = x, y = y), aes(x, y)) +
  geom_line() +
  labs(
    title = cellname,             
    x = "Time in s",              
    y = "Ca concentration in nM"  
  ) +
  theme_minimal()
p
filename = paste(cellname, "_original.png", sep = "")
ggsave(filename, p, width = 6, height = 4)


#filename = paste(cellname, "_test.png", sep = "")
peakdet_test(500, 0.025, filename)
peaks = peakdet(x, y, w=500, span=0.025)



fitted = loess(y ~ x, span = 0.025)
# Get smoothed y-values
smoothed_y = fitted$fitted
# Calculate the finite differences (approximate derivative)
dx = diff(x)
dy = diff(smoothed_y)
slope = dy / dx
x2 = x[-length(x)]
y2 = slope
fitted2 = loess(y2 ~ x2, span = 0.01)
# Get smoothed y-values
smoothed_y2 = fitted2$fitted
# Calculate the finite differences (approximate derivative)
dx2 = diff(x2)
dy2 = diff(smoothed_y2)
slope2 = dy2 / dx2
x3 = x2[-length(x2)]
y3 = slope2
fitted3 = loess(y3 ~ x3, span = 0.01)
x4 = fitted3$x
y4 = fitted3$fitted
# Calculation of possible start points
x = x4
y = y4
filename = paste(cellname, "_startpoints_test.png", sep = "")
peakdet_test(250, 0.05, filename)
points = peakdet(x, y, w=250, span=0.05)
# Look for corresponding start point for each peak 
points_x <- points$x 
peaks_x <- peaks$x
# Finding next lower values
start <- numeric(length(peaks_x))

for (i in 1:length(peaks_x)) {
  lower_values <- points_x[points_x < peaks_x[i]]
  if (length(lower_values) > 0) {
    start[i] <- max(lower_values)
  } else {
    start[i] <- NA  # No lower value found
  }
}

#reset x and y
y = original_y
x = original_x



#results table
res = data.frame(Time = x, Signal = y)
colnames(res)<-c("Time","Signal")
results_peaks = res %>%
  filter(Time %in% c(peaks$x))
results_peaks$Start_x = start
results_start = res %>%
  filter(Time %in% c(start))
colnames(results_start) = c("Start_x", "Start_y")
results = merge(results_peaks, results_start, by = "Start_x", all.x = TRUE)
colnames(results) = c("Start_x", "Peak_x", "Peak_y", "Start_y")
results = na.omit(results)
results



#find endpoints
peak_locations = results$Peak_x
peak_values = results$Peak_y
peak_baseline = results$Start_y
# Calculate the threshold as 90% of each peak's value
thresholds = peak_baseline + (peak_values - peak_baseline)*0.1
# Initialize vectors to store the x, y, and duration values
new_x = numeric(length(peak_locations))
new_y = numeric(length(peak_locations))
durations = numeric(length(peak_locations))
# Loop through each peak
for (i in 1:length(peak_locations)) {
  peak_loc = peak_locations[i]
  peak_val = peak_values[i]
  # Find the index in the x vector that corresponds to the peak location
  peak_index = which(x == peak_loc)
  # Find the index of the first point after the peak that drops below the threshold
  found_point = FALSE  # Flag to indicate if point is found
  for (j in (peak_index + 1):(length(x))) {
    if (y[j] < thresholds[i] && !found_point) {
      new_x[i] = x[j]
      new_y[i] = y[j]
      durations[i] = x[j] - peak_loc
      found_point = TRUE  # Mark the point as found
    }
  }
}
# Add new columns to the res_filtered table
results$End_x <- new_x
results$End_y <- new_y



#plot points
filename = paste(cellname, "_endpoints.png", sep = "")
png(filename, width = 8, height = 5, units = "in", res = 600)
dev.off()
dev.off()
plot(x,y,type="l",ylab="",xlab="",lwd=0.8)
lines(x = fitted$x, y = fitted$fitted, type = "l", col = "#2278AA", lwd = 2)
points(x = results$Peak_x, y = results$Peak_y , type="p", col="#EF2B2B", pch = 16)
points(x = results$Start_x, y = results$Start_y, type="p", col="#EF2B2B", pch = 4)
points(x = results$End_x, y = results$End_y , type="p", col="#EF2B2B", pch = 4)



#calculate parameters
results$Duration = results$End_x - results$Start_x
results$Amplitude = results$Peak_y - results$Start_y
# Calculate areas
areas = numeric(nrow(results))
curve_x <- original_x  # Replace with your x-values
curve_y <- original_y # Replace with your y-values
# Loop through each row of the data frame
for (i in 1:nrow(results)) {
  x = c(results$Start_x[i], results$End_x[i])  # x-coordinates of the points
  y = c(results$Start_y[i], results$End_y[i])  # y-coordinates of the points
  # Calculate the slope (m) and y-intercept (b)
  m = (y[2] - y[1]) / (x[2] - x[1])  # Slope
  b = y[1] - m * x[1]  # Y-intercept
  # Limit the curve values to the range defined by the two points
  curve_x_range = curve_x[curve_x >= x[1] & curve_x <= x[2]]
  curve_y_range = curve_y[curve_x >= x[1] & curve_x <= x[2]]
  # Calculate the limited linear line y values within the range
  limited_line_y_range = m * curve_x_range + b
  # Calculate the area using the trapezoidal rule
  areas[i] = sum(abs(curve_y_range - limited_line_y_range) * (curve_x_range[2] - curve_x_range[1]))
}
# Add the calculated areas to the data frame
results$Area = areas
filename = paste(cellname, "_final.csv", sep = "")
write.csv(results,filename)







y = original_y
x = original_x


#filter results
results_filtered = results %>%
  filter(Amplitude > 50)
results_filtered = results_filtered %>%
  filter(Start_x > 5)



##plot and save results
# Plot
filename = paste(cellname, "_final.png", sep = "")
png(filename, width = 8, height = 5, units = "in", res = 600)
plot(x,y,type="l",ylab="",xlab="",lwd=0.8) 
#lines(x = fitted$x, y = fitted$fitted, type = "l", col = "#F68989", lwd = 2)
lines(x = fitted$x, y = fitted$fitted, type = "l", col = "#2278AA", lwd = 3)
points(x = results_filtered$Peak_x, y = results_filtered$Peak_y , type="p", col="#BE0E0E", pch = 19)
points(x = results_filtered$Start_x, y = results_filtered$Start_y, type="p", col="#BE0E0E", pch = 8)
points(x = results_filtered$End_x, y = results_filtered$End_y , type="p", col="#BE0E0E", pch = 8)
dev.off()
# Save csv 
filename = paste(cellname, "_final.csv", sep = "")
write.csv(results_filtered, filename)


