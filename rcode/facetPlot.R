# facet_grid
# http://sape.inf.usi.ch/quick-reference/ggplot2/facet
# facet_grid_plot-------
d=expand.grid(obs=0:10, benchmark=c('antlr', 'bloat', 'chart', 'eclipse', 'fop', 'hsqldb', 'jython', 'luindex', 'lusearch', 'pmd', 'xalan'), gc=c('CopyMS', 'GenCopy', 'GenImmix', 'GenMS', 'Immix'), opt=c('on', 'off'), heapSize=seq(from=1.5, to=4, by=0.5))
d$time = rexp(nrow(d), 0.01)+1000
d$time = d$time + abs(d$heapSize-3)*100
d$time[d$opt=='on'] = d$time[d$opt=='on']-200

d$time[d$opt=='on' & d$benchmark=='bloat'] = d$time[d$opt=='on' & d$benchmark=='bloat'] + 190
d$time[d$opt=='on' & d$benchmark=='pmd' & d$gc=='Immix'] = d$time[d$opt=='on' & d$benchmark=='pmd' & d$gc=='Immix'] + 600
# facet_warp_plot-----
facet_wrap()
# Now let's use faceting to break the results down by benchmark:
 ggplot() +
 facet_wrap(~benchmark) +
 geom_boxplot(data=d, mapping=aes(x=opt, y=time, color=opt))
# facet_wrap(~benchmark) creates a separate panel for each benchmark. The panels
# are wrapped into multiple rows on a grid. Wrapping the panels is especially
# useful when we have a factor with a larger number of levels (such as
# benchmarks, which has 11 levels); without wrapping, the plot can become overly
# wide (or the individual panels overly narrow).
#The above plot shows that the "bloat" benchmark does not seem to benefit much
#from our optimization (its run time does not decrease much when we enable the
#optimization). The plot also shows that the optimized run time of the "pmd"
#benchmark is dispersed more (the box is taller).
# facet_grid_plot-----------------
# And now let's facet by two variables: in addition to benchmark (horizontally),
# we also use gc (vertically):
ggplot() +
facet_grid(gc~benchmark) +
geom_boxplot(data=d, mapping=aes(x=opt, y=time, color=opt))
#facet_grid(gc~benchmark) produces a separate panel for each gc-benchmark
#combination and places the panels in a grid with one row for each gc and one
#column for each benchmark. The above plot shows the reason for the wider
#dispersion of the time for the "pmd" benchmark: for most gcs, switching on the
#optimization reduces run time, however, when running on the "Immix" collector,
#the optimization actually degrades performance (increases run time).
 

 
 
 