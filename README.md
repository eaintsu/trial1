# trial1
> overall_sample_avg <- mean(plate_data$abs[!str_detect(plate_data$sample, "^Standard")])
> 
> # 2. Get all standards and calculate how far they are from your samples
> standard_selection <- plate_data %>%
+     filter(str_detect(sample, "^Standard")) %>%
+     group_by(sample) %>%
+     summarise(mean_abs = mean(abs)) %>%
+     inner_join(read_excel("/Users/eaintsu/Desktop/trial.xlsx", sheet = "Standards"), by = "sample") %>%
+     # Calculate distance from your sample average
+     mutate(dist = abs(mean_abs - overall_sample_avg)) %>%
+     # Pick the 4 standards closest to your samples
+     arrange(dist) %>%
+     slice_head(n = 4) %>%
+     arrange(std_con)
                                                                                                           
> target_model <- lm(mean_abs ~ std_con, data = standard_selection)
> t_slope     <- coef(target_model)[2]
> t_intercept <- coef(target_model)[1]
> 
> # 4. Calculate Concentrations
> final_results <- plate_data %>%
+     filter(!str_detect(sample, "^Standard")) %>%
+     group_by(sample) %>%
+     summarise(mean_sample_abs = mean(abs)) %>%
+     mutate(
+         calc_conc = (mean_sample_abs - t_intercept) / t_slope,
+         # Clean up any tiny negatives
+         calc_conc = ifelse(calc_conc < 0, 0, calc_conc)
+     )
> 
> # 5. Diagnostic Check: Which standards were used?
> print("Standards used for this curve:")
[1] "Standards used for this curve:"
> print(standard_selection)
# A tibble: 4 × 4
  sample    mean_abs std_con  dist
  <chr>        <dbl>   <dbl> <dbl>
1 Standard1    0.107       0 0.391
2 Standard2    1.19        1 0.697
3 Standard3    1.65        2 1.16 
4 Standard5    1.87       10 1.37 
> print(final_results)
# A tibble: 13 × 3
   sample mean_sample_abs calc_conc
   <chr>            <dbl>     <dbl>
 1 267B             0.314         0
 2 267L             0.496         0
 3 267R             0.406         0
 4 267RR            0.737         0
 5 270B             0.53          0
 6 270L             0.637         0
 7 270N             0.482         0
 8 270R             0.685         0
 9 274B             0.34          0
10 274L             0.683         0
11 274N             0.214         0
12 275B             0.594         0
13 275N             0.36          0
> library(readxl)
> library(tidyverse)
> library(reshape2)
> 
> # --- 1. SETTINGS ---
> file_path <- "/Users/eaintsu/Desktop/trial.xlsx"
> plate_cols <- c("alpha", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
> 
> # --- 2. LOAD & FORMAT ---
> readings <- read_excel(file_path, sheet = "Data")[, 1:13]
New names:                                                                                                 
• `` -> `...1`
• `` -> `...14`
> template <- read_excel(file_path, sheet = "Template")[, 1:13]
New names:                                                                                                 
• `` -> `...1`
> standards_key <- read_excel(file_path, sheet = "Standards")
> dilutions <- read_excel(file_path, sheet = "Dilutions")                                                  
>                                                                                                          
> colnames(readings) <- plate_cols
> colnames(template) <- plate_cols
> 
> # Melt and Merge
> readings_long <- melt(readings, id.vars = "alpha", variable.name = "col") %>% unite("well", alpha, col, sep="")
> template_long <- melt(template, id.vars = "alpha", variable.name = "col", na.rm = TRUE) %>% unite("well", alpha, col, sep="")
> plate_data <- inner_join(template_long, readings_long, by = "well") %>% setNames(c("well", "sample", "abs"))
> 
> # Calculate Means
> mean_stats <- plate_data %>% 
+     group_by(sample) %>% 
+     summarise(mean_abs = mean(abs), .groups = 'drop')
> 
> # --- 3. PLATEAU CUT-OFF & MODEL ---
> std_data <- mean_stats %>%
+     filter(str_detect(sample, "^Standard")) %>%
+     inner_join(standards_key, by = "sample") %>%
+     arrange(std_con) %>%
+     mutate(diff_abs = c(NA, diff(mean_abs)),
+            diff_con = c(NA, diff(std_con)),
+            local_slope = diff_abs / diff_con)
> 
> # Find first point where slope drops below 15% of initial segment
> initial_slope <- std_data$local_slope[2]
> linear_region <- std_data %>%
+     filter(local_slope >= (initial_slope * 0.15) | is.na(local_slope))
> 
> # Final Linear Model
> fit <- lm(mean_abs ~ std_con, data = linear_region)
> slope <- coef(fit)[2]
> intercept <- coef(fit)[1]
> 
> # --- 4. CALCULATE UNKNOWNS ---
> final_results <- mean_stats %>%
+     filter(!str_detect(sample, "^Standard")) %>%
+     mutate(calc_conc = (mean_abs - intercept) / slope,
+            calc_conc = ifelse(calc_conc < 0, 0, calc_conc)) %>%
+     left_join(dilutions, by = "sample") %>%
+     mutate(df = replace_na(df, 1),
+            final_conc = calc_conc * df)
> 
> # Output
> print(paste("R2 of Linear Region:", round(summary(fit)$r.squared, 4)))
[1] "R2 of Linear Region: 0.9477"
> print(final_results %>% select(sample, mean_abs, final_conc))
