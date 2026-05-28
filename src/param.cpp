/* param.cpp
Copyright 2009-2012 Trevor Bedford <t.bedford@ed.ac.uk>
Member function definitions for Parameter class
Parameter values are all global
*/

/*
This file is part of PACT.

PACT is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PACT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with PACT.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <Rcpp.h>

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <cstdlib>
using std::atof;

#include <stdexcept>
using std::runtime_error;
using std::out_of_range;

#include "param.h"

Parameters::Parameters(const std::string& paramFile) {

	// default parameter values
	burnin = false;
	push_times_back = false;
	reduce_tips = false;
	renew_trunk = false;
	prune_to_trunk = false;
	prune_to_time = false;
	prune_to_label = false;
	prune_to_tips = false;
	remove_tips = false;
	pad_migration_events = false;
	collapse_labels = false;
	trim_ends = false;
	section_tree = false;
	time_slice = false;
	rotate = false;
	accumulate = false;
	add_tail = false;

	print_tree = false;
	print_circular_tree = false;
	print_all_trees = false;

	summary_tmrca = false;
	summary_length = false;
	summary_root_proportions = false;
	summary_proportions = false;
	summary_coal_rates = false;
	summary_mig_rates = false;
	summary_sub_rates = false;
	summary_diversity = false;
	summary_fst = false;
	summary_tajima_d = false;
	summary_diffusion_coefficient = false;
	summary_drift_rate = false;
	summary_persistence = false;

	tips_time_to_trunk = false;
	x_loc_history = false;
	y_loc_history = false;
	coord_history = false;

	skyline_tmrca = false;
	skyline_length = false;
	skyline_proportions = false;
	skyline_coal_rates = false;
	skyline_mig_rates = false;
	skyline_pro_history_from_tips = false;
	skyline_diversity = false;
	skyline_fst = false;
	skyline_tajima_d = false;
	skyline_timetofix = false;
	skyline_xmean = false;
	skyline_ymean = false;
	skyline_xdrift = false;
	skyline_ratemean = false;
	skyline_xtrunkdiff = false;
	skyline_locsample = false;
	skyline_locgrid = false;
	skyline_drift_rate_from_tips = false;

	ordering = false;

	pairs_diversity = false;

	string paramString;
	ifstream paramFileStream(paramFile.c_str());
	if (paramFileStream.is_open()) {

		Rcpp::Rcout << "Reading parameters from " << paramFile << std::endl;
		Rcpp::Rcout << std::endl;

		while (! paramFileStream.eof() ) {
			getline(paramFileStream, paramString);
			importLine(paramString);
		}
	}
	else {
		throw runtime_error("parameter file not found: " + paramFile);
	}

}

/* Reads a string and attempts to extract parameters from it */
void Parameters::importLine(string line) {

	// READING LINE STRING
	string pstring = "";				// fill with a-z or _
	string vstring = "";				// fill with 0-9 or . or - or A-Z
	vector<double> values;				// convert vstring to double and push here
	vector<string> svalues;

	for (string::iterator is = line.begin(); is != line.end(); is++) {

		if (*is == '#')		// ignore rest of line after comment
			break;
		else if (*is >= 'a' && *is <= 'z' && vstring.size() == 0)
    		pstring += *is;
		else if ( (*is >= '0' && *is <= '9') || (*is >= 'a' && *is <= 'z') || (*is >= 'A' && *is <= 'Z') || *is == '.' || *is == '-' || *is == '/' || *is == '_' || *is == '|')
    		vstring += *is;
    	else if (vstring.size() > 0) {
    		values.push_back( atof(vstring.c_str()) );
    		svalues.push_back(vstring);
    		vstring = "";
    	}

	}

	if (vstring.size() > 0) {
		values.push_back( atof(vstring.c_str()) );
	}
	if (vstring.size() > 0) {
		svalues.push_back(vstring);
	}


	// SETTING PARAMETERS

	if (pstring == "burnin") {
		if (values.size() == 1) {
			burnin = true;
			burnin_values = values;
		}
	}

	if (pstring == "pushtimesback") {
		if (values.size() == 1 || values.size() == 2) {
			push_times_back = true;
			push_times_back_values = values;
		}
	}

	if (pstring == "reducetips") {
		if (values.size() == 1) {
			reduce_tips = true;
			reduce_tips_values = values;
		}
	}

	if (pstring == "renewtrunk") {
		if (values.size() == 1) {
			renew_trunk = true;
			renew_trunk_values = values;
		}
	}

	if (pstring == "prunetotrunk") {
		prune_to_trunk = true;
	}

	if (pstring == "prunetolabel") {
		if (values.size() == 1) {
			prune_to_label = true;
			prune_to_label_values = svalues;
		}
	}

	if (pstring == "prunetotips") {
		prune_to_tips = true;
		prune_to_tips_values = svalues;
	}

	if (pstring == "removetips") {
		remove_tips = true;
		remove_tips_values = svalues;
	}

	if (pstring == "prunetotime") {
		if (values.size() == 2) {
			prune_to_time = true;
			prune_to_time_values = values;
		}
	}

	if (pstring == "padmigrationevents") {
		pad_migration_events = true;
	}

	if (pstring == "collapselabels") {
		collapse_labels = true;
	}

	if (pstring == "trimends") {
		if (values.size() == 2) {
			trim_ends = true;
			trim_ends_values = values;
		}
	}

	if (pstring == "sectiontree") {
		if (values.size() == 3) {
			section_tree = true;
			section_tree_values = values;
		}
	}

	if (pstring == "timeslice") {
		if (values.size() == 1) {
			time_slice = true;
			time_slice_values = values;
		}
	}

	if (pstring == "rotate") {
		if (values.size() == 1) {
			rotate = true;
			rotate_values = values;
		}
	}

	if (pstring == "accumulate") {
		accumulate = true;
	}

	if (pstring == "addtail") {
		if (values.size() == 1) {
			add_tail = true;
			add_tail_values = values;
		}
	}

	if (pstring == "ordering") {
		ordering = true;
		ordering_values = svalues;
	}

	if (pstring == "printruletree" || pstring == "printtree") {
		print_tree = true;
	}

	if (pstring == "printcirculartree") {
		print_circular_tree = true;
	}

	if (pstring == "printalltrees") {
		print_all_trees = true;
	}

	if (pstring == "summarytmrca") { summary_tmrca = true; }
	if (pstring == "summarylength") { summary_length = true; }
	if (pstring == "summaryrootproportions") { summary_root_proportions = true; }
	if (pstring == "summaryproportions") { summary_proportions = true; }
	if (pstring == "summarycoalrates") { summary_coal_rates = true; }
	if (pstring == "summarymigrates") { summary_mig_rates = true; }
	if (pstring == "summarysubrates") { summary_sub_rates = true; }
	if (pstring == "summarydiversity") { summary_diversity = true; }
	if (pstring == "summaryfst") { summary_fst = true; }
	if (pstring == "summarytajimad") { summary_tajima_d = true; }
	if (pstring == "summarydiffusioncoefficient") { summary_diffusion_coefficient = true; }
	if (pstring == "summarydriftrate") { summary_drift_rate = true; }
	if (pstring == "summarypersistence") { summary_persistence = true; }

	if (pstring == "tipstimetotrunk") { tips_time_to_trunk = true; }

	if (pstring == "tipsxlochistory") {
		if (values.size() == 3) {
			x_loc_history = true;
			x_loc_history_values = values;
		}
	}

	if (pstring == "tipsylochistory") {
		if (values.size() == 3) {
			y_loc_history = true;
			y_loc_history_values = values;
		}
	}

	if (pstring == "coordhistory") {
		if (values.size() == 3) {
			x_loc_history = true;
			x_loc_history_values = values;
		}
	}

	if (pstring == "skylinesettings") {
		if (values.size() == 3) {
			skyline_values = values;
		}
	}

	if (pstring == "skylinetmrca") { skyline_tmrca = true; }
	if (pstring == "skylinelength") { skyline_length = true; }
	if (pstring == "skylineproportions") { skyline_proportions = true; }
	if (pstring == "skylinecoalrates") { skyline_coal_rates = true; }
	if (pstring == "skylinemigrates") { skyline_mig_rates = true; }
	if (pstring == "skylineprohistoryfromtips") { skyline_pro_history_from_tips = true; }
	if (pstring == "skylinediversity") { skyline_diversity = true; }
	if (pstring == "skylinefst") { skyline_fst = true; }
	if (pstring == "skylinetajimad") { skyline_tajima_d = true; }
	if (pstring == "skylinetimetofix") { skyline_timetofix = true; }
	if (pstring == "skylinexmean") { skyline_xmean = true; }
	if (pstring == "skylineymean") { skyline_ymean = true; }
	if (pstring == "skylinexdrift") { skyline_xdrift = true; }
	if (pstring == "skylineratemean") { skyline_ratemean = true; }
	if (pstring == "skylinextrunkdiff") { skyline_xtrunkdiff = true; }
	if (pstring == "skylinelocsample") { skyline_locsample = true; }
	if (pstring == "skylinelocgrid") { skyline_locgrid = true; }
	if (pstring == "skylinedriftratefromtips") { skyline_drift_rate_from_tips = true; }

	if (pstring == "pairsdiversity") {
		if (values.size() == 1) {
			pairs_diversity = true;
			pairs_diversity_values = values;
		}
	}

}

/* prints parameters */
void Parameters::print() {

	// GENERAL
	if ( general() ) {

		Rcpp::Rcout << "General:" << std::endl;

		if (burnin) {
			Rcpp::Rcout << "burnin " << burnin_values[0] << std::endl;
		}

		Rcpp::Rcout << std::endl;

	}

	// TREE MANIPULATION
	if ( manip() ) {

		Rcpp::Rcout << "Tree manipulation:" << std::endl;

		if (push_times_back) {
			Rcpp::Rcout << "push times back ";
			for (int i = 0; i < (int)push_times_back_values.size(); i++) {
				Rcpp::Rcout << push_times_back_values[i] << " ";
			}
			Rcpp::Rcout << std::endl;
		}

		if (reduce_tips) {
			Rcpp::Rcout << "reduce tips " << reduce_tips_values[0] << std::endl;
		}

		if (renew_trunk) {
			Rcpp::Rcout << "renew trunk " << renew_trunk_values[0] << std::endl;
		}

		if (trim_ends) {
			Rcpp::Rcout << "trim ends " << trim_ends_values[0] << " " << trim_ends_values[1] << std::endl;
		}

		if (section_tree) {
			Rcpp::Rcout << "section tree " << section_tree_values[0] << " " << section_tree_values[1] << " " << section_tree_values[2] << std::endl;
		}

		if (time_slice) {
			Rcpp::Rcout << "time slice " << time_slice_values[0] << std::endl;
		}

		if (prune_to_label) {
			Rcpp::Rcout << "prune to label " << prune_to_label_values[0] << std::endl;
		}

		if (prune_to_tips) {
			Rcpp::Rcout << "prune to tips:" << std::endl;
			for (int i = 0; i < (int)prune_to_tips_values.size(); i++) {
				Rcpp::Rcout << prune_to_tips_values[i] << " ";
			}
			Rcpp::Rcout << std::endl;
		}

		if (remove_tips) {
			Rcpp::Rcout << "remove tips:" << std::endl;
			for (int i = 0; i < (int)remove_tips_values.size(); i++) {
				Rcpp::Rcout << remove_tips_values[i] << " ";
			}
			Rcpp::Rcout << std::endl;
		}

		if (prune_to_trunk) {
			Rcpp::Rcout << "prune to trunk" << std::endl;
		}

		if (prune_to_time) {
			Rcpp::Rcout << "prune to time " << prune_to_time_values[0] << " " << prune_to_time_values[1] << std::endl;
		}

		if (pad_migration_events) {
			Rcpp::Rcout << "pad migration events" << std::endl;
		}

		if (collapse_labels) {
			Rcpp::Rcout << "collapse labels" << std::endl;
		}

		if (rotate) {
			Rcpp::Rcout << "rotate " << rotate_values[0] << std::endl;
		}

		if (accumulate) {
			Rcpp::Rcout << "accumulate" << std::endl;
		}

		if (add_tail) {
			Rcpp::Rcout << "add tail " << add_tail_values[0] << std::endl;
		}

		if (ordering) {
			Rcpp::Rcout << "Tip ordering:" << std::endl;
			for (int i = 0; i < (int)ordering_values.size(); i++) {
				Rcpp::Rcout << ordering_values[i] << " ";
			}
			Rcpp::Rcout << std::endl;
		}

		Rcpp::Rcout << std::endl;

	}

	// TREE STRUCTURE
	if ( printtree() ) {
		Rcpp::Rcout << "Tree structure:" << std::endl;
		if (print_tree) { Rcpp::Rcout << "print tree" << std::endl; }
		if (print_circular_tree) { Rcpp::Rcout << "print circular tree" << std::endl; }
		if (print_all_trees) { Rcpp::Rcout << "print all trees" << std::endl; }
		Rcpp::Rcout << std::endl;
	}

	// SUMMARY STATISTICS
	if ( summary() ) {
		Rcpp::Rcout << "Summary statistics:" << std::endl;
		if (summary_tmrca) { Rcpp::Rcout << "tmrca" << std::endl; }
		if (summary_length) { Rcpp::Rcout << "length" << std::endl; }
		if (summary_root_proportions) { Rcpp::Rcout << "root proportions" << std::endl; }
		if (summary_proportions) { Rcpp::Rcout << "proportions" << std::endl; }
		if (summary_coal_rates) { Rcpp::Rcout << "coal rates" << std::endl; }
		if (summary_mig_rates) { Rcpp::Rcout << "mig rates" << std::endl; }
		if (summary_sub_rates) { Rcpp::Rcout << "sub rates" << std::endl; }
		if (summary_diversity) { Rcpp::Rcout << "diversity" << std::endl; }
		if (summary_fst) { Rcpp::Rcout << "fst" << std::endl; }
		if (summary_tajima_d) { Rcpp::Rcout << "tajima d" << std::endl; }
		if (summary_persistence) { Rcpp::Rcout << "persistence" << std::endl; }
		Rcpp::Rcout << std::endl;
	}

	// TIP STATISTICS
	if ( tips() ) {
		Rcpp::Rcout << "Tip statistics:" << std::endl;
		if (tips_time_to_trunk) { Rcpp::Rcout << "time to trunk" << std::endl; }
		if (x_loc_history) {
			Rcpp::Rcout << "x loc history " << x_loc_history_values[0] << " " << x_loc_history_values[1] << " " << x_loc_history_values[2] << std::endl;
		}
		if (y_loc_history) {
			Rcpp::Rcout << "y loc history " << y_loc_history_values[0] << " " << y_loc_history_values[1] << " " << y_loc_history_values[2] << std::endl;
		}
		Rcpp::Rcout << std::endl;
	}

	// SKYLINE STATISTICS
	if ( skyline() ) {
		Rcpp::Rcout << "Skyline statistics:";
		Rcpp::Rcout << " " << skyline_values[0];
		Rcpp::Rcout << " " << skyline_values[1];
		Rcpp::Rcout << " " << skyline_values[2] << std::endl;
		if (skyline_tmrca) { Rcpp::Rcout << "tmrca" << std::endl; }
		if (skyline_length) { Rcpp::Rcout << "length" << std::endl; }
		if (skyline_proportions) { Rcpp::Rcout << "proportions" << std::endl; }
		if (skyline_coal_rates) { Rcpp::Rcout << "coal rates" << std::endl; }
		if (skyline_mig_rates) { Rcpp::Rcout << "mig rates" << std::endl; }
		if (skyline_pro_history_from_tips) { Rcpp::Rcout << "pro history from tips" << std::endl; }
		if (skyline_diversity) { Rcpp::Rcout << "diversity" << std::endl; }
		if (skyline_fst) { Rcpp::Rcout << "fst" << std::endl; }
		if (skyline_tajima_d) { Rcpp::Rcout << "tajima d" << std::endl; }
		if (skyline_timetofix) { Rcpp::Rcout << "time to fix" << std::endl; }
		if (skyline_xmean) { Rcpp::Rcout << "x mean" << std::endl; }
		if (skyline_ymean) { Rcpp::Rcout << "y mean" << std::endl; }
		if (skyline_xdrift) { Rcpp::Rcout << "x drift" << std::endl; }
		if (skyline_ratemean) { Rcpp::Rcout << "rate mean" << std::endl; }
		if (skyline_xtrunkdiff) { Rcpp::Rcout << "x trunk diff" << std::endl; }
		if (skyline_locsample) { Rcpp::Rcout << "loc sample" << std::endl; }
		if (skyline_locgrid) { Rcpp::Rcout << "loc grid" << std::endl; }
		if (skyline_drift_rate_from_tips) { Rcpp::Rcout << "drift rate from tips" << std::endl; }
		Rcpp::Rcout << std::endl;
	}

	// PAIR STATISTICS
	if ( pairs() ) {
		Rcpp::Rcout << "Pair statistics:" << std::endl;
		if (pairs_diversity) { Rcpp::Rcout << "pairwise diversity" << std::endl; }
	}

}

bool Parameters::general() {
	return burnin;
}

bool Parameters::manip() {
	return (push_times_back || reduce_tips || renew_trunk || prune_to_trunk || prune_to_label ||
	        prune_to_tips || remove_tips || pad_migration_events || collapse_labels || trim_ends ||
	        section_tree || time_slice || rotate || accumulate || add_tail || ordering);
}

bool Parameters::printtree() {
	return (print_tree || print_circular_tree || print_all_trees);
}

bool Parameters::summary() {
	return (summary_tmrca || summary_length || summary_root_proportions || summary_proportions ||
	        summary_coal_rates || summary_mig_rates || summary_sub_rates || summary_diversity ||
	        summary_fst || summary_tajima_d || summary_diffusion_coefficient || summary_persistence);
}

bool Parameters::tips() {
	return (tips_time_to_trunk || x_loc_history || y_loc_history || coord_history);
}

bool Parameters::skyline() {
	return ( skyline_values.size() == 3 &&
	        (skyline_tmrca || skyline_length || skyline_proportions || skyline_coal_rates ||
	         skyline_mig_rates || skyline_pro_history_from_tips || skyline_diversity ||
	         skyline_fst || skyline_tajima_d || skyline_timetofix || skyline_xmean ||
	         skyline_ymean || skyline_xdrift || skyline_ratemean || skyline_xtrunkdiff ||
	         skyline_locsample || skyline_locgrid || skyline_drift_rate_from_tips) );
}

bool Parameters::pairs() {
	return pairs_diversity;
}
