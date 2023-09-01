### Initial Quality Control (light version)
(called via initialQualityControl_light.m)

Performing the following checks:
- for fields that have to be present
- Number of metabolites, reactions and genes
- Number of subsystems & Number of reactions per subsystem
- Values of lower and upper bounds
- Reactions in objective function
- flux consistency
- Number and ratio of consistent reactions
- Number of blocked/closed reactions
- Maximal growth rate
- Minimal mode for enabling growth and respective number of reactions per subsystem
- Minimal medium for enabling growth

### Initial Quality Control
(called via callQualityControl.m; needs fastbox toolbox)

Performing the following checks:
- for fields that have to be present
- if identifiers are unique
- for empty rows in the S matrix
- if rules and GPR rules are consistent
- for duplicated reactions
- for metabolites without formulas
- for dead-ends
- if subsystem are given
- if reactions are assigned to multiple subsystems

Furthermore:
- reports flux consistency and blocked/closed reactions
- simplify reversibility
- reports stats on metabolites, reactions and genes
- reports minimal bound_span (ub-lb) to identify narrow bounded reactions
- reports number of optimized reactions
- reports minimal mode for enabling growth and respective number of reactions per subsystem
- reports respective minimal medium for enabling growth
