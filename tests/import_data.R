# import_data.R
# create a test dataset for import data feature
d = data.frame(
  y = c(0.1, 0.2, 0.5, 0.7, 0.8, 0.8),
  x = c(1, 2, 3, 4, 5, 6)
)
write.csv(d, "tests/import_test_data/import_test.csv", row.names = F)
