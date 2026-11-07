## code to prepare `expr_toy` dataset goes here

# 生成模拟数据
set.seed(123)
genes <- c("NSUN2", "TG", "TSHR", "PAX8", "SLC5A5", "NKX2-1", paste0("Gene", 1:94))
expr_toy <- matrix(rnorm(300), nrow = 100, ncol = 3, dimnames = list(genes, paste0("Sample", 1:3)))
expr_toy["Gene1", ] <- expr_toy["NSUN2", ] * 0.8 + rnorm(3, sd = 0.2)
expr_toy["Gene2", ] <- -expr_toy["TG", ] * 0.8 + rnorm(3, sd = 0.2)

# 保存进 R 包
usethis::use_data(expr_toy, overwrite = TRUE)
