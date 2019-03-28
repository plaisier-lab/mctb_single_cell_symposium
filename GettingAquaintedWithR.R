##########################################################
## MCTB Single Cell Symposium: GettingAquaintedWithR.R  ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
##                                                      ## 
## @Course:  MCTB Single Cell Symposium                 ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @github: https://github.com/plaisier-lab/            ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

# Set repositories
setRepositories()

# What is my working directory
getwd()

# Set your working directory
#setwd() do this through the GUI

# R as a calculator
1 + 1

# Create a vector from 1 to 10
c(1,2,3,4,5,6,7,8,9,10)

# An easier way to create a vector from 1 to 10
1:10

# Make a vector of 10 1's
rep(1,10)

# What is rep? Putting a question mark in front of a function pulls up it's help page.
?rep

# Assigning an object to a variable
a1 <- rep(1,10)

# Take a look at the object
a1

# Again
b1 <- 1:10
b1

# Concatenating (combining) vectors
c(a1, b1)

# Math with vectors (Note:  arithmetic operations act element-wise on vectors.)
a1+b1
a1+2*b1
a1/b1

# Accessing elements of a vector
b1
b1[4]
b1[2:4]

# Subset based upon a threshold
b1<7
b1[b1<7]

# Sorting of vectors
c1 = c(8,5,6,3,2,4,9,10,1)
sort(c1)

# Make a matrix
m1 = matrix(nrow = 10, ncol = 3, data = c(1:10,10:1,1:10))
m1

# Accessing rows from a matrix
m1[1,]
m1[10,]

# Access columns from a matrix
m1[,1]
m1[,2]
m1[,3]

# Sort matrix by column 2
order(m1[,2])
m1[order(m1[,2]),]

# What vignettes are available
vignette()

# Opening a vignette
vignette('regression-examples')

# Play with vignette code
vig = vignette('regression-examples')
edit(vig)


## Simple for Loop
# P-values from our analysis
p.values = c(0.1, 0.05, 0.003, 0.4, 0.9)
# A vector to store the negative log p-values
neglog10.p.values = 1:5
# Transform the p-values
for(p in 1:length(p.values)) {
  neglog10.p.values[p] = -log10(p.values[p])
}


## Simple while loop
# Numbers from our analysis
v1 = c(21, 22, 53, 74, 85, 96, 97, 58, 49, 30, 85)
# Iterator
i = 1
# Look for first instance of 85
while(v1[i] != 85) {
  i = i + 1 
}
# Print out where we found it
print(paste('v1[',i,'] = 85',sep=''))


## Simple functions - Which is faster for or while?
# Multiply numbers together for all elements of N using for loop
fn1 <- function(N) {
  for(i in as.numeric(1:N)) {
    y <- i*i
  }
}

# Multiply numbers together for all elements of N using while loop
fn2 <- function(N) {
  i = 1
  while(i <= N) {
    y <- i*i
    i <- i + 1
  }
}

# How long do they take
system.time(fn1(60000))
system.time(fn2(60000))


## T-test
t.test(sleep[1:10,'extra'],sleep[11:20,'extra'])
# Or using formula
t.test(extra ~ group, data = sleep)

## Correlation
cor.test(c(1:10), c(11:20))
cor.test(c(1:10), -c(11:20))

# What else is built-in
library(help="stats")

