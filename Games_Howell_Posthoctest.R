# Lade benötigte Pakete
library(ggpubr)
library(agricolae)
library(car)
library(stats)
install.packages("jmv")
library(jmv)


# Überprüfe die Varianz-Homogenität
levene_test_result <- leveneTest(Read_Length ~ Name, data = combined_cores_ts)
print(levene_test_result)

# Überprüfe die Normalverteilung
# Funktion zum Erstellen eines QQ-Plots für eine Stichprobe
create_qqplot <- function(group_data, group_name) {
  sample_data <- sample(group_data$Read_Length, size = min(100, length(group_data$Read_Length)), replace = FALSE)
  ggplot(data.frame(Read_Length = sample_data), aes(sample = Read_Length)) +
    geom_qq() +
    geom_qq_line() +
    ggtitle(group_name)
}

# Liste der QQ-Plots für jede Gruppe
qqplots_list <- by(combined_cores_ts, combined_cores_ts$Name, create_qqplot, group_name = names(combined_cores_ts$Name))

# Kombiniere die QQ-Plots
gridExtra::grid.arrange(grobs = qqplots_list, ncol = 1)

# Beispiel: Erstelle einen Datensatz (ersetze dies mit deinen eigenen Daten)
set.seed(42)  # Damit die Zufallszahlen reproduzierbar sind
combined_cores_ts <- rnorm(5000)

# Führe den Shapiro-Wilk-Test durch
shapiro_test <- shapiro.test(combined_cores_ts)

# Zeige das Ergebnis des Tests an
print(shapiro_test)

# Führe den Welch's-Test durch
welchs_test <- oneway.test(Read_Length ~ Name, data = combined_cores_ts)

# Zeige das Ergebnis des Tests an
print(welchs_test)


# trying with jmv package
result <- anovaOneW(formula = Read_Length ~ Name,
                    data = combined_cores_ts,
                    welchs = TRUE,
                    norm = TRUE,
                    eqv = FALSE,
                    phMethod = 'gamesHowell')

posthoc_results <- as.data.frame(result$postHoc)

write.csv(posthoc_results, file = "posthoc_results.csv", row.names = FALSE)


