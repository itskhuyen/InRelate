
library(testthat)
library(Rsolnp)
source("InRelate.R")

#This test reads the example.indivq file and checks if the readIndivq function correctly extracts the subpopulation information (Pop column) from it.


# Read the example.indivq file
indivq_path <- "../example.indivq"
example_indivq <- read.table(indivq_path, header = FALSE, sep = "\t")

test_that("readIndivq works correctly with example.indivq", {
  result <- readIndivq(indivq_path)
  expect_equal(result, example_indivq)
  subpops <- getSubpop(result)
  expect_equal(subpops,3)
})


# Read the example.pkla file
pkla_path <- "../examplepkla.txt"
example_pkla <- read.table(pkla_path, header = FALSE, sep = "\t")

test_that("readPkla works correctly with example.pkla", {
  result <- readPkla(pkla_path)
  expect_equal(result, example_pkla)
})

# Read the example.str file
str_path <- "../example.str"
example_str <- read.table(str_path, header = FALSE, sep = "\t",skip=1)

test_that("readStr works correctly with example.str", {
  result <- readStr(str_path)
  expect_equal(result, example_str)
})


#contigency test 1
test_that("check that datatset is a dataframe", {
  expected_data_type <- "data.frame"
  result <- readStr(str_path)
  expect_s3_class(result, expected_data_type)
})


#contigency test 1: the sum of the values after the fifth column adds up to 1
test_sum_of_example_indivq <- function(example_indivq){
  result <- sum(example_indivq[, 6:ncol(example_indivq)]) == 1
  expected <- TRUE
  if (!identical(result, expected)){
    stop("Test Failed")
  }
}

#contigency test 2
test_that("check that datatset is a dataframe", {
  expected_data_type <- "data.frame"
  result <- readStr(indivq_path)
  expect_s3_class(result, expected_data_type)
})


#contigency test 1
test_that("check that datatset is a dataframe", {
  expected_data_type <- "data.frame"
  result <- readStr(pkla_path)
  expect_s3_class(result, expected_data_type)
})

#contigency test 2/ check that across columns after the second for one locus adds up to 1
sum_of_unique_values <- function(example_pkla, k) {
  unique_value <- unique(example_pkla[, "V1"])
  result <- data.frame(V1 = unique_value)
  for (j in 2:k) {
    colmn_sum <- numeric()
    for (i in 1:length(unique_value)) {
      colmn_sum[i] <- sum(example_pkla[example_pkla[, "V1"] == unique_value[i], j])
    }
    result[, j] <- colmn_sum
  }
  return(result)
}

#test_that("sum of unique values in V1 is equal to 1", {
#  expected <- sum_of_unique_values(example_pkla)
#  result <-
#})

#test for getting unique individuals from the example_str
#unit test
test_that("get unique indiv", {
  result <- getIndiv(example_str)
  expected <- as.numeric(unique(example_str$V1))
  expect_equal(result,expected)
})

#test for getting K subpopulation from indinvq file
#unit test
#test_that("get k subpop", {
#  result <- getSubpop(example_str)
#  expected <- example_str[,6:ncol(example_str)]
#  expect_equal(result,expected)
#})

#check that the file is a dataframe
test_that("check dataframe", {
  result <- getLoci(example_str)
  expect_true(is.data.frame(example_str))
})

#test get unique alleles from example_pkla
#unit test
test_that("get unique alleles", {
  result <- getAlleles(example_pkla)
  expected <- aggregate(x = example_pkla$V2, by = list(example_pkla$V1), FUN = function(x) length(unique(x)))
  expect_equal(result,expected)
})

#check that the file is a dataframe
test_that("check dataframe", {
  result <- getAlleles(example_pkla)
  expect_true(is.data.frame(example_pkla))
})

#test for getting loci from example_str
#unit test
#test_that("get loci loci", {
#  result <- getLoci(example_str)
#  expected <- example_str[,3:ncol(example_str)]
#  expect_equal(result,expected)
#})
 
#check that the file is a dataframe
test_that("check dataframe", {
  result <- getLoci(example_str)
  expect_true(is.data.frame(example_str))
})

#test functions for Fst function from Pegas
# test is using Jaquar dataset from R package
#test to check for diploid input
#test_that("check for non-diploid input", {
#  x <- jaguar
#  expect_error(Fst(x), "Fst() requires diploid data input")
#})

#check for population assignment
#test_that("check for population assignment", {
#  x <- jaguar
#  expect_equal(Fst(x, pop = as.factor(4)), Fst(x, pop = as.factor(4)))
#})

#test function for computenumpairs
#indiv function use in computenumpairs
getIndiv <- function(example_str){
  individual <- as.numeric(unique(example_str$V1));
  return(individual)
}


#check that getIndiv function is reporting the correct individuals
test_that("get unique indiv", {
  result <- getIndiv(example_str)
  expected <- as.numeric(unique(example_str$V1))
  expect_equal(result,expected)
})







