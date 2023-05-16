library(ggplot2)
ggplot(data = mtcars, aes(x = wt, y = mpg)) +
  annotate("text", x = 2.5, y = 25, 
           label = "bgroup('(',atop(x,y),')')", 
           parse = TRUE) +
  annotate("text", x = 3.5, y = 25, 
           label = "group('(',atop(x,y),')')", 
           parse = TRUE)


ggsave(filename = "bgroup rendering issue\\test_1.png"
       plot = last_plot())


ragg::agg_png(filename = "bgroup rendering issue\\test_5.png", )
plot(0, xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression(bgroup('(',atop(x,y),')')))
dev.off()
dev.off()
