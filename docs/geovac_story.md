# GeoVac in plain language

*A short account of what this project is, for a reader with no background in the field. Every statement here is made
precisely, with its evidence, in the papers; this page only tells the story. Nothing in the project has been
externally peer-reviewed.*

## It started with counting

Put points on rings so that they are evenly spaced, and let the rings grow. The numbers of points that fit go
1, 3, 5, 7. Those are also the numbers of orbitals in the s, p, d and f families that organise the periodic table.
The match is in the counts. Calling the third ring "p" is a label we choose, not something the rings tell us.

## There is a sphere inside the hydrogen atom

In 1935 Vladimir Fock showed that the hydrogen atom, looked at in the right way, is a wave on the surface of a
four-dimensional ball. The atom's energy levels are the ways a wave can wrap around that surface, and the labels
physicists give to electron states (n, l, m) are the labels of those wrappings.

The counting above gives exactly that list of labels. So the project's starting picture, a lattice of labelled
points, and the ordinary equation of quantum mechanics are two descriptions of the same thing. One limit is known and
recorded: the lattice has the right labels and the right energies, but the links between its points do not turn into
the sphere as the lattice grows.

## Leaving the sphere has a price

On the lattice everything is a whole number or a fraction. Numbers like pi and logarithms appear only when you leave
it: when you measure, add a second atom, or add up over all of space. The project's working claim is that each way of
leaving brings in its own kind of number, and that you can say which.

Mathematicians call the numbers that come out of such sums *periods*, and sort them by how complicated they are. The
project keeps a list, called the atlas, with one row for each way of leaving the sphere: what breaks, what survives,
and what kind of number it costs.

## The ladder

Sorted this way, the building blocks of chemistry sit on a ladder.

- **One atom:** fractions, and fractions times a square root.
- **Two atoms:** logarithms and their relatives.
- **Three atoms:** numbers that come from curves shaped like the surface of a doughnut.
- **More atoms, bending, squashing, outside fields:** in every case checked so far, nothing higher. Changing the
  geometry alone seems to stop at the doughnut.
- **Linking more electrons together:** this keeps climbing, one step for each electron cloud in the chain.

The higher steps are the same mathematical objects that particle physicists meet when they calculate collisions.
That is the project's main finding: the mathematics of molecules and the mathematics of particle physics share a
ladder, and the project has mapped where the molecular pieces sit on it, with proofs for some steps and computer
checks for others.

## What it is not

It is not a new theory of nature. The sphere and the usual equations say the same thing.

It is not a replacement for the chemistry software people use today. The project does build compact descriptions of
molecules for quantum computers, and they are sparse and quick to make, but they are rough. In the one careful
head-to-head on accuracy (water, with a tuned basis), the advantage over standard methods was real but modest, and
its size depends on which standard method you compare with.

## How it is done

The project is run by one independent researcher working with AI models. That makes care more important, not
less. Predictions are written down before the calculation is run. Each result is attacked by a fresh reviewer, and
gets a test that has been shown to fail when the answer is wrong. Dead ends are kept on the record with the reason
they failed. Most of what the map shows was known in pieces to specialists, and the papers say so; what is new is
drawing it in one place for chemistry and checking it this carefully.

## Where to go next

- The project as a longer story: `papers/synthesis/geovac_field_guide.pdf`
- The ladder, for readers with a mathematics or physics background: `papers/synthesis/period_ladder.pdf`
- Every headline claim with its evidence: `docs/claims_register.md`
- The list of ways out of the sphere: `docs/atlas.md`
