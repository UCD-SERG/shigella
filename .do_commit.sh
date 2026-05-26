#!/bin/bash
git -c user.name="Claude" -c user.email="claude[bot]@users.noreply.github.com" \
  commit --author="claude[bot] <claude[bot]@users.noreply.github.com>" \
  -m "fix(sampler): increase max_treedepth default from 12 to 15

50% treedepth saturation on n=5 phase0 run (500/1000 transitions hit
the ceiling) is the dominant driver of ESS_bulk=16 and R-hat=1.766.
The LKJ Cholesky geometry requires longer HMC trajectories than
max_treedepth=12 allows; raising to 15 gives 8x more trajectory
length budget before truncation.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
