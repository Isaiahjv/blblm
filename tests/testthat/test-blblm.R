test_that("type test", {

  data1 = data.frame(x = rnorm(100,20,5), z = rnorm(100,80,9), y = rpois(100,4))
  data2 = data.frame(x = rnorm(100,20,5), z = rnorm(100,80,9), y = rpois(100,4))

  m = blblm(y~., data = data1, model = 'poisson', m = 3, B = 100, core = 1, para = TRUE)
  coef = blbcoef(m)

  expect_s3_class(m, "blblm")

  #testing near equality with glm function. randomness pulls it away
  expect_equal(coef(glm(y~.,data = data1, family = 'poisson')), coef.blblm(m),tolerance=.6)

  #testing near equality of prediction
  expect_equal(predict(glm(y~.,data = data1, family = 'poisson'),data2),predict.blblm(m,data2),tolerance=.6)



})
