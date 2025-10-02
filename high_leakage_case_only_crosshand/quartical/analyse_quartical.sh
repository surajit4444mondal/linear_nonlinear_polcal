echo "Calibrating S1"


goquartical calibrate_unpolarized_source.yaml

echo "Calibrating S2"

goquartical calibrate_polarized_Q1_U5_V0.yaml

echo "Calibrating S3"

goquartical calibrate_polarized_Q5_U1_V0.yaml

echo "Calibrating S4"

goquartical polarized_source_I50_Q4_U15_V0.yaml

echo "Calibrating S5"

goquartical polarized_source_I50_Q15_neg_U5.yaml
