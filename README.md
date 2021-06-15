# To build
`rivet-build SSHiggs.cc`

# To run
`rivet --analysis=SSHiggs --pwd /eos/user/s/sdysch/SSHiggs/outputs/schan_ZGamma_ll_eemumu_MH++_default/tag_1_pythia8_events.hepmc.gz --nominal-weight Weight_MERGING=0.000 --skip-weights`

# To make an html plotbook
`rivet-mkhtml -c SSHiggs.plot -o plots SSHiggs.yoda:"Title:Output"`
