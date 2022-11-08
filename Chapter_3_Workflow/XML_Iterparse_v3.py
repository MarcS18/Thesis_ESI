# -*- coding: utf-8 -*-
from lxml import etree
import csv

"""
XML parser for extracting chemical reaction data using organometallic catalysts from a .xml file extracted
from the Reaxys chemical database. This example is designed explicitly for the Ullmann-Goldberg reaction (copper) with
accompanying filters.
"""


def main():

    # Import xml file and input name for output file
    xmlfile = input("Enter input file name")
    outfile = input("Enter output file name")
    xml_parse(xmlfile, outfile)


def xml_parse(xmlfile, outfile):

    # Define catalyst filter terms
    cat_metal = ['copper', 'Cu', 'Copper', 'cuprous']
    # Define general reagents for filter
    genrgts = ['potassium', 'sodium', 'caesium', 'cesium', 'aluminum', 'aluminium', 'hydrogenchloride', 'lithium',
               'calcium', 'magnesium', 'barium', 'hydrogen chloride', 'silicate', 'crown', 'ammonium', 'imidazolium',
               'iron', 'palladium', 'Raney', 'Ni', 'phosphazene', 'silica', 'poly', 'alumina', 'fluoride', 'sieve',
               'cyanide', 'methoxide', 'NaOCH3', 'triethylamine', 'sulfuric', 'trifluoroacetic acid', 'nitrogen',
               'carbonate', 'hydrogen',  'thionyl', 'silver', 'jones', 'nitric', 'permanganate', 'oxygen',
               'selenium', 'mercury', 'phosphorus', 'choramine', 'boron', 'citric acid', 'ethylamine', 'tin',
               'sulfonic acid', 'iodine',  'trichlorophosphate', 'sulfonate', 'oxalyl', 'zinc', 'transferase',
               'NADPH', 'sulphurous acid', 'ammomium', 'NaH', 'nickel', 'K2CO3', 'hypophosphorous acid', 'silane',
               'p-methoxybenzyl chloride', 'Grignard', 'K3PO4', 'nitrite',  'titanium', 'air', 'bismuth', 'arsenic',
               'hydrazine hydrate', 'manganese', 'methyl iodide', 'chromate', 'carbon dioxide', 'rhodium', 'ammonim',
               'iridium', 'Fe', 'HATU', 'formic acid', 'K2O3', 'nitrous acid', 'cobalt', 'Cs2O3', 'Na,K-tartrate']
    # Define solvents for filter
    gen_solv = ['ethyl acetate', 'dimethyl sulfoxide', 'DMSO', 'propan-1-ol', 'Ethyl tert-butyl ether',
                'pentan-1-ol', 'i-Amyl alcohol', 'N,N-dimethyl-formamide', 'ethanol', 'toluene', 'butan-1-ol', 'benzene',
                 'hexan-1-ol', 'naphthalene', 'tert-butyl alcohol', 'xylene', ' dimethyl acetamide',
                'nitrobenzene', '1-methyl-pyrrolidin-2-one', 'nitromethane', '1,2-dichlorobenzene', 'NMP', 'acetonitrile',
                'tetrahydrofuran', 'diethyl ether', 'dichloromethane', 'isopropyl alcohol', 'chloroform', '1,4-dioxane',
                'Isopropyl acetate', 'water', 'hexadecane', 'ethyl methyl ether', '1,3,5-trimethyl-benzene',
                'dodecane', 'DMF', 'diphenylether', 'Diphenylmethane', 'dimethylamine', 'dihexyl ether', 'dibutyl ether',
                'decalin', 'decane', '5,5-dimethyl-1,3-cyclohexadiene', '2-methoxy-ethanol', '2-ethoxy-ethanol'
                , '1,2-dimethoxyethane', 'hexachloroethane', 'Isobutyronitrile', 'paraffin', 'phenanthrene', 'propyl methanoate'
                , 'formamide', 'Petroleum ether', 'ethylene glycol', 'cyclohexanone']

    # Remove any bold, italic, subscript, superscript, highlights from xml file
    with open(xmlfile) as f:
        xml_str = f.read()

    xml_str = xml_str.replace('<hi>', '')
    xml_str = xml_str.replace('</hi>', '')
    xml_str = xml_str.replace('<sub>', '')
    xml_str = xml_str.replace('</sub>', '')
    xml_str = xml_str.replace('<SUB>', '')
    xml_str = xml_str.replace('</SUB>', '')
    xml_str = xml_str.replace('<sup>', '')
    xml_str = xml_str.replace('</sup>', '')
    xml_str = xml_str.replace('<SUP>', '')
    xml_str = xml_str.replace('</SUP>', '')
    xml_str = xml_str.replace('<i>', '')
    xml_str = xml_str.replace('</i>', '')
    xml_str = xml_str.replace('<br/>', '')

    with open(xmlfile, 'w') as f:
        f.write(xml_str)

    # Assign lists
    data = []
    reactants = []
    rct1 = []
    rct2 = []
    reagents = []
    all_reagents = []
    all_catalysts = []
    product = []
    prep = []
    Yield = []
    stage = []
    steps = []
    catalyst = []
    ligand = []
    solvent = []
    time = []
    temperature = []
    ref = []

    # Column titles
    Titles = ['Index', 'Reactant 1', 'Reactant 2', 'Product', 'Product MW', 'Procedure', 'Yield (%)', 'Steps', 'Reagent(s)', 'Ligand(s)',
              'Catalyst(s)', 'Solvent(s)', 'Time (h)', 'Temperature (°C)', 'All Reagents', 'All Catalysts', 'Reference Type'
        , 'Reference']

    # Write column titles to file
    with open(outfile, 'w+', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(Titles)

    # Parse XML file
    for event, element in etree.iterparse(xmlfile, events=("start", "end")):
        # Reaction index
        if element.tag == "reaction" and event == "start":
            reaction = element.attrib['index']
        # Append reactants to list
        elif element.tag == "RX.RCT" and event == "start":
            reactant = element.text
            reactants.append(reactant)
        # Assign reaction product
        elif element.tag == "RX.PRO" and event == "start":
            product = element.text
        # Assign product molecular weight
        elif element.tag == "RX.MAXPMW" and event == "start":
            PMW = element.text
        # Assign reaction procedure
        elif element.tag == "RXD.TXT" and event == "start":
            prep = element.text
        # Assign reaction yield
        elif element.tag == "RXD.NYD" and event == "start":
            Yield = element.text
        # Append number of reaction steps
        elif element.tag == "RXD.STG" and event == "start":
            stage.append(element.text)
        # Append reagents to all_reagents list and pass each through catalyst, solvent and reagent filter
        # else assign as ligand
        elif element.tag == "RXD.RGT" and event == "start":
            reagent = element.text
            all_reagents.append(reagent)
            if reagent is None:
                pass
            else:
                if any(r in reagent for r in cat_metal):
                    catalyst = reagent
                elif any(r in reagent for r in genrgts):
                    reagents.append(reagent)
                elif any(r in reagent for r in gen_solv):
                    solvent.append(reagent)
                else:
                    ligand = reagent
        # Append catalyst to all_catalysts list and pass through catalyst, reagent and solvent filter
        # else assign as ligand
        elif element.tag == "RXD.CAT" and event == "start":
            cat = element.text
            all_catalysts.append(cat)
            if cat is None:
                pass
            else:
                if any(r in cat for r in cat_metal):
                    catalyst = cat
                elif any(r in cat for r in genrgts):
                    reagents.append(cat)
                elif any(r in cat for r in gen_solv):
                    solvent.append(cat)
                else:
                    ligand = cat
        # Append solvent to list
        elif element.tag == "RXD.SOL" and event == "start":
            solvent.append(element.text)
        # Append reaction time to list
        elif element.tag == "RXD.TIM" and event == "start":
            time.append(element.text)
        # Append temperature to list
        elif element.tag == "RXD.T" and event == "start":
            temperature.append(element.text)
        # Assign citation type
        elif element.tag == "CIT.DT" and event == "start":
            ref_type = element.text
        # Append patent author to reference list
        elif element.tag == "CIT.PA" and event == "start":
            ref.append(element.text)
        # Append patent year to reference list
        elif element.tag == "CIT.PREPY" and event == "start":
            ref.append(element.text)
        # Append patent pages to reference list
        elif element.tag == "CIT.PAGES" and event == "start":
            ref.append(element.text)
        # Append page number to reference list
        elif element.tag == "CIT.PN" and event == "start":
            ref.append(element.text)
        # Append author to reference list
        elif element.tag == "CIT.AU" and event == "start":
            ref.append(element.text)
        # Append journal to reference list
        elif element.tag == "CIT.JTS" and event == "start":
            ref.append(element.text)
        # Append volume number to reference list
        elif element.tag == "CIT.VL" and event == "start":
            ref.append(element.text)
        # Append issue number to reference list
        elif element.tag == "CIT.NB" and event == "start":
            ref.append(element.text)
        # Append publication year to reference list
        elif element.tag == "CIT.PY" and event == "start":
            ref.append(element.text)
        # Append page number to reference list
        elif element.tag == "CIT.PAG" and event == "start":
            ref.append(element.text)
        # Append DOI to reference list
        elif element.tag == "CIT.DOI" and event == "start":
            ref.append(element.text)
        # End of <RXD> entry write to file
        elif element.tag == "RXD" and event == "end":
            # Separate reactants
            rct1 = reactants[0]
            try:
                rct2 = reactants[1]
            except:
                rct2 = "None"
            # Assign number of reaction steps
            if len(stage) == 0:
                steps = '1'
            else:
                steps = stage[-1]
            # Convert all lists to strings with | separator
            rgt =" | ".join([str(x) for x in reagents])
            all_rgts = " | ".join([str(x) for x in all_reagents])
            all_cats = " | ".join([str(x) for x in all_catalysts])
            solvs = " | ".join([str(x) for x in solvent])
            times = " | ".join([str(x) for x in time])
            temperatures = " | ".join([str(x) for x in temperature])
            ref2 = ", ".join([str(x) for x in ref])

            # Write variables to list
            data = [reaction, rct1, rct2, product, PMW, prep, Yield, steps, rgt, ligand, catalyst, solvs, times, temperatures, all_rgts, all_cats, ref_type, ref2]
            # If no value assign blank
            data = ['None' if i == '' else i for i in data]
            # Write data to file
            with open(outfile, 'a+', encoding='utf-8', newline='') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerow(data)
            # Clear all variables and lists in <RXD>
            rct1 = ""
            rct2 = ""
            prep = ""
            Yield = ""
            steps = ""
            rgt = ""
            all_rgts = ""
            all_cats = ""
            catalyst = ""
            ligand = ""
            ref.clear()
            stage.clear()
            reagents.clear()
            all_reagents.clear()
            all_catalysts.clear()
            solvent.clear()
            time.clear()
            temperature.clear()
            data.clear()
        # At end of reaction clear reactants, products and molecular weight
        elif element.tag == "reaction" and event == "end":
            reactants.clear()
            product = ""
            PMW = ""
            element.clear()
        # Else remove from memory
        else:
            pass
            element.clear()

    # Correct/Remove characters in output file
    with open(outfile) as f:
        csv_str = f.read()

    csv_str = csv_str.replace('(l)', '(I)')
    csv_str = csv_str.replace('`', "'")
    csv_str = csv_str.replace('[]', "None")
    csv_str = csv_str.replace(';', "")

    with open(outfile, 'w') as f:
        f.write(csv_str)


if __name__ == '__main__':
    main()
