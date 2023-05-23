# Example to make publication quality tables

library(tables)
library(tablesgg)
iris2_tab <- tabular(Species*Heading()*value*Format(digits=2)*(mean + sd) ~
                       Heading("Flower part")*flower_part*Heading()*direction,
                     data=iris2)
plot(textTable(iris2_tab))
