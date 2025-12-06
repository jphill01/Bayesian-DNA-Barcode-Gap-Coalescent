library(rstudioapi)

# Define your folder paths
folder1 <- selectDirectory()
folder2 <- selectDirectory()
folder3 <- selectDirectory()


# Function to extract species names from subfolder names
get_species_2 <- function(folder) {
  dirs <- list.dirs(folder, full.names = FALSE, recursive = FALSE)
  dirs[dirs != ""]  # remove empty entries
}

# Collect species from each folder
species1 <- get_species_2(folder1)
species2 <- get_species_2(folder2)
species3 <- get_species_2(folder3)

# All unique species
all_species <- unique(c(species1, species2, species3))

# Partition species into mutually exclusive groups
only_f1   <- setdiff(species1, union(species2, species3))
only_f2   <- setdiff(species2, union(species1, species3))
only_f3   <- setdiff(species3, union(species1, species2))

f1_f2     <- setdiff(intersect(species1, species2), species3)
f1_f3     <- setdiff(intersect(species1, species3), species2)
f2_f3     <- setdiff(intersect(species2, species3), species1)

all_three <- Reduce(intersect, list(species1, species2, species3))

# Print results
cat("Total number of unique species:", length(all_species), "\n\n")

cat("1. Species only in folder1:\n")
print(only_f1)

cat("\n2. Species only in folder2:\n")
print(only_f2)

cat("\n3. Species only in folder3:\n")
print(only_f3)

cat("\n4. Species in folder1 & folder2 only:\n")
print(f1_f2)

cat("\n5. Species in folder1 & folder3 only:\n")
print(f1_f3)

cat("\n6. Species in folder2 & folder3 only:\n")
print(f2_f3)

cat("\n7. Species in all three folders:\n")
print(all_three)

# Sanity check: all groups together should equal all_species
partitioned <- c(only_f1, only_f2, only_f3, f1_f2, f1_f3, f2_f3, all_three)
cat("\nSanity check passed? ", setequal(all_species, partitioned), "\n")

