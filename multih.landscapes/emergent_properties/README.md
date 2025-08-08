# Are food web triads simply the sum of their habitats?

The analyses developped here aim at testing whether the network properties we measured on sites combining multiple habitats are different than expected if landscape scale food webs were only the sum of their habitats.

This test is conducted by means of null models used to generate random plant-insect interaction networks in triads under the hypothesis that triad networks are simply the sum of their habitats.

To that end, these null models follow at least the following rules:

* Random triad networks (for a given site) are the result of three random networks for each habitat of that triad. These random networks are random subsamples of the corresponding monads, and should preserve either the number of interaction events or the level of sampling completeness achieved in each habitat of the triad.

* When subsampling monads to generate random sub-triad networks, interaction events observed in monads are sampled with a probability that is proportional to the number of times they have been observed during the sampling survey.

The present folder is structured with three sub-folders: one for the R scripts necessary to this analysis (`main_scripts`), one for the functions (`functions`), and one for the outputs such as tables and figures (`outputs`).
