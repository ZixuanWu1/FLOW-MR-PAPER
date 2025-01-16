df = read.csv("bmiA_on_sbpA.csv")
ggplot(df, aes(x = Estimate, y = Method)) + xlim(-0.05, 0.2) +
  geom_point() +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.2) +
  labs(y = "Method", x = "Estimate") +
  geom_vline(xintercept = 0,linetype="dotted")+
  theme_bw()
