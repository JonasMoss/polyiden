#The local correlation should have be a gaussian that matches up with the
#given box-probabilities. However, there are two solutions.
#Several possible ways to fix this could be attempted.
#One that may seem acceptable is the following:

#Example:
#
# __________
#|     |    |
#| A   |  B |
#|_____|____|tau22
#|     |    |
#| C   | D  |
#|_____|____|tau21
#|     |    |
#| E   | F  |
#|_____|____|
#tau11, tau12

# For box E, we have by the polychoric case a unique correlation attached to it
#from matching CDFs at the top-right corner.

#For box C, there are likely two solutions. However, we also have a unique solution
#from the "sister" problem of matching the CDF at C's top most right corner.
#Perhaps it makes sense to choose the solution to the "box-problem" that is
#closest to this unique CDF-matching?

#Another idea is: For each box, we can take the correlation that is the average of
#the CDF in each corner where we have a CDF-restriction. So for B, it is the prob
#P ( Z_1 > \tau_12, Z_2 > \tau_22), and this is the only restriction (the rest
#are given by the marginals). For box C we have two CDF-restrictions.
#For more categories, some boxes would have four.

#Each square is compatible with normality. However, it does not need to be the case
#that the square has restrictions compatible with the other CDF-restrictions connected to that
#box.