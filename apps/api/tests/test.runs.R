context("Testing all runs endpoints")

test_that("Calling /api/runs/ with a valid workflow id returns Status 200", {
  res <- httr::GET(
    "http://localhost:8000/api/runs/?workflow_id=1000009172",
    httr::authenticate("carya", "illinois")
  )
  expect_equal(res$status, 200)
})



test_that("Calling /api/runs/{id} with a valid run id returns Status 200", {
  res <- httr::GET(
    "http://localhost:8000/api/runs/1002042201",
    httr::authenticate("carya", "illinois")
  )
  expect_equal(res$status, 200)
})

test_that("Calling /api/runs/ with a invalid workflow id returns Status 404", {
  res <- httr::GET(
    "http://localhost:8000/api/runs/?workflow_id=1000000000",
    httr::authenticate("carya", "illinois")
  )
  expect_equal(res$status, 404)
})



test_that("Calling /api/runs/{id} with a invalid run id returns Status 404", {
  res <- httr::GET(
    "http://localhost:8000/api/runs/1000000000",
    httr::authenticate("carya", "illinois")
  )
  expect_equal(res$status, 404)
})