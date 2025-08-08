* `structural-calculations_mdt-MANOVA.R`
	* This is the code to extract diversity and structural metrics and analyse the effect of number of habitats on community diversity and structure. It also runs a MANOVA and provides post-hoc analyses.
	* Input files are stored in `../../data` in this repository:
		* `all_web_interactions_final.txt`
		* `site_level_data.txt`
		* `MDT_floral-diversity.csv`
		* `extracted_site_level_data.txt` which corresponds to all extracted metric outputs and site data generated in 1 file.
	* Outputs are gathered in `outputs`:
		* `insects_total-number_per-site.txt`
		* `insects_total-sp-richness_per-site.txt`
		* `floral-units_per-sites.txt`
		* `insects-plants_sp-evenness.txt`
		* `all-webs_interaction-evenness.txt`

* `structural_analysis_figure2.R`
	* It generates figure 2, a 6-panel figure for floral abundance, plant species richness, insect abundance, insect species richness, insect species evenness and interaction evenness for monads, dyads and triads.
	* `extracted_site_level_data.txt` is the input file for this figure.
