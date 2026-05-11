# Set working directory manually first if needed
# For Permethrin
perm <- prop.test(22, 100)
# For Deltamethrin
delt <- prop.test(98, 100)

#Check the range under the 95 percent confidence interval
perm
delt