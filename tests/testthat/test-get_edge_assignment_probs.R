test_that("get_edge_assignment_probs works", {
    iseas <- get_edge_assignment_probs(m = 8, deg.seq = c(4, 4, 4, 4), model = "IEAS")
    isa <- get_edge_assignment_probs(m = 10, deg.seq = c(8, 4, 2, 2, 2, 2), model = "ISA")

    iseas_correct <- c(
        0.05, 0.133333333333333, 0.133333333333333, 0.133333333333333,
        0.05, 0.133333333333333, 0.133333333333333, 0.05, 0.133333333333333,
        0.05
    )
    isa_correct <- c(
        0.16, 0.16, 0.08, 0.08, 0.08, 0.08, 0.04, 0.04, 0.04, 0.04,
        0.04, 0.01, 0.02, 0.02, 0.02, 0.01, 0.02, 0.02, 0.01, 0.02, 0.01
    )

    expect_equal(sort(iseas), sort(iseas_correct))
    expect_equal(sort(isa), sort(isa_correct))

    expect_error(get_edge_assignment_probs(m = 8, deg.seq = c(1, 2, 4, 4), model = "IEAS"))
    expect_error(get_edge_assignment_probs(m = 8, deg.seq = c(4, 4, 4, 4), model = "BLABLA"))
})
