


p1 <- ggplot() +
  theme_bw() +
  geom_ribbon(aes(x = ID, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`)), data = m_2$summary.random[[2]], alpha = 0.3) +
  geom_line(aes(x = ID, y = exp(mean)), data = m_2$summary.random[[2]]) +
  ylim(0,8) +
  xlab("Sea Ice Concentration (%)") +
  ggtitle("Effect of sea ice concentration on harp seal occurence") +
  geom_vline(xintercept = 15, linetype = "dashed") +
  geom_vline(xintercept = 80, linetype = "dashed")

quartz(width = 6, height = 5)
print(p1)
quartz.save("~/harp seals and sea ice.jpeg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()


p2 <- ggplot() +
  geom_ribbon(aes(x = ID, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`)), data = m_2$summary.random[[3]], alpha = 0.3) +
  geom_line(aes(x = ID, y = exp(mean)), data = m_2$summary.random[[3]]) +
  xlab("Deviation from seasonal average sea ice concentration (%)") +
  coord_cartesian(ylim = c(0, 20))

quartz(width = 10, height = 6)
gridExtra::grid.arrange(p1, p2, ncol = 2)
quartz.save("2 component ice output.jpeg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()

require(animation)
saveVideo({
  for (i in 1:20){
    p1 <- ggplot() +
      theme_bw() + ylab("") + xlab("") +
      gg(mesh, col = m_1$summary.random$s$mean[idx$s.group == i]) + #+ m_1$summary.random$`inla.group(ice_av, n = 25, method = "cut")`$mean) +
      scale_fill_viridis("", limits = c(-5, 25), breaks = seq(-5, 25, 5), na.value = "transparent") +
      geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
      coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
      ggtitle(i)
    print(p1)
  }
}, movie.name = "INLA_harps_20.mp4", interval = 0.5, ani.width = 750, ani.height = 750, other.opts = "-pix_fmt yuv420p -b:v 1080k")




m_1$summary.random$s$mean[idx$s.group == i]
m_1$summary.random$`inla.group(ice_av, n = 25, method = "cut")`



# The seasonal pattern has been lost by shifting from 4 season cyclic model to 20 season continous model
# How do I add seasonality to the continous model?
# Will running with more data help? 1000 locations rather than 500 locations?


ggplot() +
  theme_bw() + ylab("") + xlab("") +
  gg(mesh)

gg(mesh, col = m_1$summary.random$s$mean) +
  facet_wrap(~ idx$s.group)



foo <- inla.stack.data(stk, tag = "pred")
preds <- m_1$summary.linear.predictor