# Hu_Code
 2012 Hu code and upgrades (C)

Note that while I have tried to benchmark the templates in this code to the Hu+2012 paper, the results do not always perfectly reproduce that paper. In particular, the H2 scenario has very different CO, CH4, and somewhat different upper-altitude H2O behaviour. I attribute this to the 2012 paper being a moving target, and it being up in the air what version of the code was used for the different scenarios (and what the inputs were). For example, the H2 PTZ grid in what Renyu sent us originally was different from the one in the input files he provided us. The irradiation may also have changed, as may the network (both the H2 and N2 templates had the wrong NKinM parameter (smaller), indicating they were written for an older version of the code.)

In maintaining the code, when I make changes, I re-run the 8 scenarios from the Hu+2012 paper, i.e. Earth, Mars, CO2-dominated atmosphere with varying levels of outgassing orbiting sun, Earthlike planet orbiting Sun with varying atmospheric redox states. I also do the exoplanet with varying redox states orbiting Sun with the extended chemical network (including N species).
