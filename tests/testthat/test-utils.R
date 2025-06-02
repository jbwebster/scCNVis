
#These are very basic and cover very few edge cases

test_that("Can create scCNVis object", {
  obj <- createPlotObject(example.matrix, example.cells, example.granges, example.meta)
  expect_equal(class(obj) == "scCNVisObject", TRUE)
})


test_that("Invalid scCNVis object creation is detected", {
  expect_error(createPlotObject("A", example.cells, example.granges, example.meta))
  expect_error(createPlotObject(example.matrix, "A", example.granges, example.meta))
  expect_error(createPlotObject(example.matrix, example.cells, "A"))
})


test_that("AddMetaData() is functional", {
  obj <- createPlotObject(example.matrix, example.cells, example.granges, example.meta)
  obj <- addMetaData(obj, c(1:length(example.cells)), "Example1")
  expect_equal(length(colnames(obj@Meta)), 2)
  expect_error(addMetaData(obj, c(1:4000), "Example2")) #Too many values
  objx <- addMetaData(obj, "X", "Example2", rownames(obj@Matrix)[1])
  expect_equal(nrow(obj@Meta), nrow(objx@Meta))
})

