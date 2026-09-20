# Adapter repair v3

The first managed run (162bd494-205b-4ce2-abc7-8ad09f740533) failed before inference because the input adapter read the record carrier instead of its payload. Both consuming modules now use the SDK unwrap_payload helper. A real BioWorld wiring regression exercises the three nodes with an explicitly simulated CLI failure; test fixtures are never benchmark predictions. The BioSimulant dependency is pinned to 0.0.34 to match the observed managed runtime. Scientific inputs, Boltz 2.0.2, sampling controls, analysis and acceptance criteria are unchanged. The original failed attempt remains in the report.

## Adapter repair v4
Run 04061d28-b96f-4258-9346-a4f794572817 reached provenance collection but failed because the managed interpreter does not include pip. No durable scientific artifact was returned, so no prediction is counted. Package versions now use importlib.metadata; hardware inspection is bounded and nonfatal. Scientific inputs and analysis are unchanged.
