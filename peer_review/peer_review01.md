## Peer Review 01

**Reviewer:** Kai Liu

**Reviewee:** Cenying Yang

**Date:** Sep 18th, 2016


Overall, both R code and mathematical derivation for exercises 01 and 02 are
great. They are very clean and easy to follow. I actually went over Google's R
style guide, and found that you almost follow all rules. For example, do not
use underscores or hyphens in both variable and function names; function names
have initial capital letters and no dots. However, there are also some other
style requirements that you should follow:

- The maximum line length is 80 characters. It makes readers read your code
  more easily since they do not need to scroll left and right over and over
  again.

- Use `<-`, not `=`, for assignment. (Actually I did not follow this rule as well)

- Functions should contain a comments section immediately below the function
  definition line, including description of the function, a list of the
  function's arguments, and a description of the return value.
  ([Reference](https://google.github.io/styleguide/Rguide.xml#identifiers))

In addition, it would be better if you could post visualized results of those
code on GitHub.

Also, two specific comments:

- In functions for gradient descent algorithm and Newton's method, you used
  threshold to constrain the function which is great, especially when the
  algorithm converges before reaching the maximum iterations. But what if it
  does not converge after finishing maximum iterations? You have no way to
  track it... So I recommend to return both betas and log likelihood values in
  the function so that you can track the convergence of the algorithm.

- In function for stochastic gradient descent algorithm, you put initial betas
  into it. The downside here is that you cannot change them when running the
  function. I think it will be better if you could treat betas as an argument
  of the function.

I also learn one function from your code:

- function `solve`. I always use to compute the inverse of a matrix, but never
  realize that it can be used to calculate x in a function `Ax=b`.
