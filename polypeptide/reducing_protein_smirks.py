import os
from openeye import oechem
import pickle
import copy
import glob
import itertools
from chemper.chemper_utils import check_smirks_to_reference, get_typed_molecules, create_tuples_for_clusters
from making_proteins import everything_from_fastas, print_order_type_data, at_least_one_passed, ParameterSystem, by_biggest_size, by_biggest_smirks
import cmiles
import json
from parmed.modeller import ResidueTemplate
from chemper.smirksify import Reducer, print_smirks


def convert_json_and_oeb(json_file, mol_dir='./mol_files/'):
    """
    Takes a json file created during the making_proteins step.
    Finds the oeb files associated with that dictionary
    and extracts the OEMol objects.

    Returns
    -------
    mols: list of OEMols
    smirks_dict: dictionary with automatically created SMIRKS for
                 different fragment types
    clusters: atom indice clusters used to make the SMIRKS patterns
    """
    with open(json_file, 'r') as inputf:
        d = json.load(inputf)

    mol_dir = os.path.abspath(mol_dir)
    mol_files = [os.path.join(mol_dir, m) for m in d['mol_files']]
    mols = list()
    for mol_file in mol_files:
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(mol_file)
        while oechem.OEReadMolecule(ifs,mol):
            mols.append(oechem.OEMol(mol))

    return mols, d['smirks_lists'], d['clusters']


# create dictionary to store final SMIRKS
final_dict = dict()

# list to find all relevant files
file_keys = [
    ('big', ['big_smirks', 'biggest_size']),
    ('small', ['small_smirks', 'small_size'])
]

# loop over both big and small order types
for fn_label, cluster_orders in file_keys:
    # find files for that order type
    fns = glob.glob('./mol_files/allIn1_%s_99sbildn_*_1mols.json' % fn_label)
    print('='*80)
    print(' '*20,fn_label)
    print('='*80)
    for f in fns:
        frag = f.split('_')[-2]
        # if its a torsion parameter type find proper or improper
        if frag == 'torsion':
            prefix = f.split('_')[-3]
            frag = '%s_%s' % (prefix, frag)
        print('-'*80)
        print(' '*30,frag)
        print('-'*80)
        # convert json file to get SMIRKS and molecules
        mols, dsmirks, dclusters = convert_json_and_oeb(f)

        # if we haven't made a dictionary for this fragment
        # make a subdictionary
        if frag not in final_dict:
            final_dict[frag] = dict()

        for order in cluster_orders:
            final_dict[frag][order] = dict()
            d = dsmirks[order][frag]
            type_list = [(l, s) for l,s in d['type_list']]
            final_dict[frag][order]['initial'] = type_list

            print('ORIGINAL', order)
            print_smirks(type_list)

            if not d['checked']:
                # wasn't able to make SMIRKS for this order
                # and fragment type combination
                # note: this means there's an incorrect key in the dict.
                #       I forgot to switch this to output_5k instead of 10
                final_dict[frag][order]['output_10k'] = None
                continue

            # make a Reducer to generalize the SMIRKS patterns
            red = Reducer(type_list, mols, verbose=False)

            # Run 1k iterations and save and print the SMIRKS
            final_dict[frag][order]['output_1k'] = red.run(1000)
            print('REDUCED 1k ', order)
            print_smirks(final_dict[frag][order]['output_1k'])

            # Run 4k more iterations (5k total)
            # save and print the SMIRKS
            final_dict[frag][order]['output_5k'] = red.run(4000)
            print('REDUCED 5k ', order)
            print_smirks(final_dict[frag][order]['output_5k'])



# Pickle dictionary for posterity
pickle.dump(final_dict, open('./mol_files/reduced_smirks_dict_5k.p', 'wb'))


