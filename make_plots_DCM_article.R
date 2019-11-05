### FIGURE 1

# Show raw chla, raw cdom and raw doxy ==> don't show corrected chla because the FDOM correction was not introduced in the text yet
# we could take a float with the deep offset correction BUT they don't have CDOM fluorometers..

tmp <- profiles[profiles$platform == 6901866,]
tmp <- tmp[tmp$id == floor(mean(unique(tmp$id))),]
tmp_doxy <- doxy[doxy$juld == unique(tmp$juld),]

# Normalization of all data to avoid multiples X axes (good practice?)
old_chla_norm <- normalize(tmp$old_fluo)
cdom_norm <- normalize(tmp$cdom)

# Replace true data with normalized data
tmp$old_fluo <- old_chla_norm
tmp$cdom <- cdom_norm


ggplot(tmp, aes(x=old_fluo, y = depth)) + geom_path(colour = "green4") + scale_y_reverse() + theme_bw() +
  xlab(expression(Normalized~data)) + ylab("Depth (m)") + theme(text=element_text(size=12)) +
  geom_path(data = tmp, aes(x=cdom, y = depth), colour = "red") + 
  geom_path(data = tmp_doxy, aes(x = doxy_norm, y = depth), colour = "black")

### END OF FIGURE 1
