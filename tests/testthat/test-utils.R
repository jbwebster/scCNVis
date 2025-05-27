
#These are very basic and cover very few edge cases

test_that("Can create SCCNVis object", {
  obj <- createPlotObject(example.matrix, example.cells, example.granges)
  expect_equal(class(obj) == "SCCNVisObject", TRUE)
})


test_that("Invalid SCCNVis object creation is detected", {
  expect_error(createPlotObject("A", example.cells, example.granges))
  expect_error(createPlotObject(example.matrix, "A", example.granges))
  expect_error(createPlotObject(example.matrix, example.cells, "A"))
})


test_that("AddMetaData() is functional", {
  obj <- createPlotObject(example.matrix, example.cells, example.granges)
  obj <- addMetaData(obj, c(1:length(example.cells)), "Example1")
  expect_equal(length(colnames(obj@Meta)), 2)
  expect_error(addMetaData(obj, c(1:4000), "Example2")) #Too many values

})

