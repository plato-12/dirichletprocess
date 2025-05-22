# inst/dev/dependency_graph.R
library(DiagrammeR)

# Create a graph of implementation dependencies
create_dependency_graph <- function() {
  grViz("
  digraph implementation_plan {
    # Node styling
    node [shape=box, style=filled, fillcolor=lightblue, fontname=Arial, fontsize=12];

    # Core data structures
    subgraph cluster_0 {
      label = 'Core Data Structures';
      style = filled;
      color = lightgrey;

      base_classes [label='Base Classes'];
      distribution_classes [label='Distribution Classes'];
      special_classes [label='Special Case Classes'];

      base_classes -> distribution_classes;
      base_classes -> special_classes;
    }

    # Core MCMC components
    subgraph cluster_1 {
      label = 'Core MCMC Components';
      style = filled;
      color = lightgrey;

      likelihood [label='Likelihood Functions'];
      ccu [label='ClusterComponentUpdate'];
      cpu [label='ClusterParameterUpdate'];
      update_alpha [label='UpdateAlpha'];
      mh [label='Metropolis-Hastings'];

      likelihood -> ccu;
      likelihood -> cpu;
      likelihood -> mh;
      mh -> cpu;
    }

    # Utility functions
    subgraph cluster_2 {
      label = 'Utility Functions';
      style = filled;
      color = lightgrey;

      rng [label='Random Number Generation'];
      predictive [label='Predictive Distributions'];
      cluster_mgmt [label='Cluster Management'];

      rng -> predictive;
      predictive -> ccu;
      ccu -> cluster_mgmt;
    }

    # Integration
    subgraph cluster_3 {
      label = 'R Integration';
      style = filled;
      color = lightgrey;

      interface [label='R/C++ Interface'];
      memory [label='Memory Management'];

      base_classes -> interface;
      distribution_classes -> interface;
      interface -> memory;
    }

    # Testing
    subgraph cluster_4 {
      label = 'Testing & Validation';
      style = filled;
      color = lightgrey;

      unit_tests [label='Unit Tests'];
      integration_tests [label='Integration Tests'];

      ccu -> unit_tests;
      cpu -> unit_tests;
      update_alpha -> unit_tests;
      interface -> integration_tests;
    }

    # Cross-cluster dependencies
    distribution_classes -> likelihood;
    distribution_classes -> rng;
    ccu -> interface;
    cpu -> interface;
    update_alpha -> interface;
  }
  ")
}

# Save the graph
create_dependency_graph()
