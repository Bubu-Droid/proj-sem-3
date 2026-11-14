# Ordinal and Latent-Grade Modelling of Radiologist-Assessed Lung Nodule Malignancy

This was the project topic for my third semester statistics project.

- PDF: [main.pdf](main.pdf)
- TeX: [main.tex](main.tex)
- Slides: [stat-presentation-final.pdf](stat-presentation-final.pdf)

You can find the dataset used in this project from [TCIA][lidcidri].

## Project Abstract

This project models radiologists' malignancy ratings of lung nodules in the
LIDC-IDRI database, using geometric features computed from their outlines. We
compare a proportional odds model of the median rating with a latent-grade model
fitted by the EM algorithm. Size is the dominant predictor, but the proportional
odds assumption is rejected. On held-out patients, the models predict exact grades
poorly and the EM model does not improve on always predicting grade $`3`$. However, a
triage rule flags $`91.3\%`$ of high-grade nodules while flagging $`47.0\%`$ of all
nodules.

All models are fitted to the trimmed file
`nodules_train_trimmed.csv` rather than the entire training data: we
removed individual ratings at least two grades from their nodule's median ($`212`$
of the $`4647`$ training ratings, $`4.6\%`$, across $`204`$ nodules) and recomputed each
median. One nodule lost all of its ratings and is excluded from the naive models.
The test set uses all ratings. This rule was chosen after an initial fit gave
implausible results. It improves the interpretability of the latent classes, but
it reduces measured disagreement, so the estimated confusion matrix overstates
reader reliability, and predictive fit on the untrimmed ratings is worse.

[lidcidri]: https://www.cancerimagingarchive.net/collection/lidc-idri/
