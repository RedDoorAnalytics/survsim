# survsim

## Simulating time-to-event data from parametric distributions, custom distributions, competing risk models and general multi-state models

For an introductory seminar, see [Simulating time-to-event data from parametric distributions, custom distributions, competing risk models and general multi-state models](https://www.youtube.com/watch?v=hmeE0qPOjP8).

## Installation

The current stable version of `survsim` can be installed from the SSC archive by typing:

```{stata}
ssc install survsim
```

To install directly from this GitHub repository, use:

```{stata}
net install survsim, from("https://raw.githubusercontent.com/RedDoorAnalytics/survsim/main/")
```

# References

> Bender R, Augustin T and Blettner M. Generating survival times to simulate Cox proportional hazards models. *Statistics in Medicine* 2005;24:1713-1723.

> Beyersmann J, Latouche A, Buchholz A and Schumacher M. Simulating competing risks data in survival analysis. *Statistics in Medicine* 2009;28:956-971.
    
> Crowther MJ and Lambert PC. Simulating complex survival data. *The Stata Journal* 2012;12(4):674-687.

> Crowther MJ and Lambert PC. Simulating biologically plausible complex survival data. *Statistics in Medicine* 2013;32(23):4118-4134.

> Crowther MJ. Simulating time-to-event data from parametric distributions, custom distributions, competings risk models and general multi-state models. Pre-print 2020.

> Jann, B. 2005. moremata: Stata module (Mata) to provide various functions. Available from http://ideas.repec.org/c/boc/bocode/s455001.html.
