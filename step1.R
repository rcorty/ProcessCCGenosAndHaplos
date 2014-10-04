library(data.table)
library(ggplot2)

file.names <- list.files(path = 'CC_genos_and_haplos_from_web', full.names = FALSE)
file.paths <- list.files(path = 'CC_genos_and_haplos_from_web', full.names = TRUE)
num.files <- length(file.names)

for (file.num in seq(from = 1, to = num.files)) {
  
  file.name <- file.names[file.num]
  file.path <- file.paths[file.num]
  strain.name <- strsplit(file.name, split = '_')[[1]][1]
  
  print(strain.name)
  
  one.strain.dt <- data.table(strain.name = strain.name,
                              read.csv(file = file.path, stringsAsFactors = FALSE))
  setnames(one.strain.dt, c('strain.name', 'chr', 'pos', 'genotype', 'haplotype'))
  
  if (file.num == 1) {
    all.strain.dt <- one.strain.dt
  } else {
    all.strain.dt <- rbind(all.strain.dt, one.strain.dt)
  }
}

# make sure things make sense from a strain perspective
unique(all.strain.dt$strain.name)
all.strain.dt[, .N, by = strain.name]
markers.by.strain.plot <- ggplot(data = all.strain.dt, aes(x = strain.name)) + 
  geom_bar(binwidth = 1, color = 'blue', fill = 'darkgray') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = 'CC_markers_by_strain.png', plot = markers.by.strain.plot, width = 12, height = 8)


# make sure things make sense from a chr perspective
unique(all.strain.dt$chr)
barplot(table(all.strain.dt$chr))
all.strain.dt[chr == 'X', chr := '20']
all.strain.dt[, chr := as.numeric(chr)]
unique(all.strain.dt$chr)
markers.by.chr.plot <- ggplot(data = all.strain.dt, aes(x = chr)) + 
  geom_bar(binwidth = 1, color = 'blue', fill = 'darkgray') + 
  xlim(c(1,20))
ggsave(filename = 'CC_markers_by_chr.png', plot = markers.by.chr.plot, width = 12, height = 8)

# make sure position column makes sense
marker.density.plot <- ggplot(data = all.strain.dt, aes(x = pos)) + 
  geom_histogram() + 
  facet_wrap(~chr) + 
  ggtitle('CC marker density across genome') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = 'CC_marker_density.png', plot = marker.density.plot, width = 12, height = 8)

# make sure things make sense wrt haplotype
all.strain.dt[, .N, by = haplotype]
haplotypes.prevalence.plot <- ggplot(data = all.strain.dt, aes(x = haplotype)) + 
  geom_bar(binwidth = 1, color = 'blue', fill = 'darkgray')
ggsave(filename = 'CC_haplotype_prevalence.png', plot = haplotypes.prevalence.plot, width = 12, height = 9)

# replace founder strains with their letters
all.strain.dt[haplotype == 'A/J', haplotype := 'A']
all.strain.dt[haplotype == 'C57BL/6J', haplotype := 'B']
all.strain.dt[haplotype == '129S1/SvImJ', haplotype := 'C']
all.strain.dt[haplotype == 'NOD/LtJ', haplotype := 'D']
all.strain.dt[haplotype == 'NZO/HILtJ', haplotype := 'E']
all.strain.dt[haplotype == 'CAST/EiJ', haplotype := 'F']
all.strain.dt[haplotype == 'PWK/PhJ', haplotype := 'G']
all.strain.dt[haplotype == 'WSB/EiJ', haplotype := 'H']
unique(all.strain.dt$haplotype)

# save the genotype and haplotype info
saveRDS(object = all.strain.dt, file = 'CC_genos_and_known_haplos.RDS')
