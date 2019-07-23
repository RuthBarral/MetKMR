test_annotation <- data.frame(

)

test_that("description", {
    # Reorder the original data randomnly
    data_mod <- testthat_data[sample(1:nrow(testthat_data)), ]

    # Create first object with original data
    instance_org <- new("MetRKAT",
                        data = testthat_data,
                        wsize = 1000,
                        wmethod = "location",
                        distmethod = c("euclidean", "manhattan"),
                        gap = 400,
                        max.na = 0.5,
                        annotation = sim_annotation_df)


})

test_that("invalid window size", {
    expect_error(new("MetRKAT",
                        data = testthat_data,
                        wsize = -2,
                        wmethod = "location",
                        distmethod = c("euclidean", "manhattan"),
                        gap = 400,
                        max.na = 0.5,
                        annotation = testthat_annot))
})

test_that("invalid gap size", {
    expect_error(new("MetRKAT",
                     data = testthat_data,
                     wsize = 100,
                     wmethod = "location",
                     distmethod = c("euclidean", "manhattan"),
                     gap = 400,
                     max.na = 0.5,
                     annotation = testthat_annot))
})

test_that("invalid distmethod size", {
    expect_error(new("MetRKAT",
                     data = testthat_data,
                     wsize = 100,
                     wmethod = "location",
                     distmethod = "invaliddist",
                     gap = 400,
                     max.na = 0.5,
                     annotation = testthat_annot))
})

test_that("description", {

})
