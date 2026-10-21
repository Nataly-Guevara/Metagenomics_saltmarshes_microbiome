x <- read.csv("assembly_process_pairwise_with_metadata.csv",
              check.names = FALSE)

make_summary <- function(data, group_var) {
  tab <- table(data[[group_var]], data$Process)

  out <- as.data.frame(tab)
  names(out) <- c("Comparison", "Process", "Count")

  totals <- aggregate(Count ~ Comparison, out, sum)
  names(totals)[2] <- "Total"

  out <- merge(out, totals, by = "Comparison")
  out$Percentage <- round(100 * out$Count / out$Total, 2)

  out[order(out$Comparison, -out$Percentage), ]
}

same_site <- make_summary(x, "Same_site")
same_treatment <- make_summary(x, "Same_treatment")
same_compartment <- make_summary(x, "Same_compartment")

site_pairs <- make_summary(x, "Site_comparison")
treatment_pairs <- make_summary(x, "Treatment_comparison")
compartment_pairs <- make_summary(x, "Compartment_comparison")

write.csv(same_site,
          "assembly_percent_same_vs_different_site.csv",
          row.names = FALSE)

write.csv(same_treatment,
          "assembly_percent_same_vs_different_treatment.csv",
          row.names = FALSE)

write.csv(same_compartment,
          "assembly_percent_same_vs_different_compartment.csv",
          row.names = FALSE)

write.csv(site_pairs,
          "assembly_percent_by_site_pair.csv",
          row.names = FALSE)

write.csv(treatment_pairs,
          "assembly_percent_by_treatment_pair.csv",
          row.names = FALSE)

write.csv(compartment_pairs,
          "assembly_percent_by_compartment_pair.csv",
          row.names = FALSE)

cat("\nSame versus different sites:\n")
print(same_site)

cat("\nSame versus different treatments:\n")
print(same_treatment)

cat("\nSame versus different compartments:\n")
print(same_compartment)

cat("\nSite-pair percentages:\n")
print(site_pairs)
