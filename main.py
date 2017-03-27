from __future__ import print_function
import warnings
import sys
import os

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.File import UndoHandle
from Bio import BiopythonWarning
from Bio.Alphabet import generic_protein
from SCOPData import protein_letters_3to1

global mmb_path

mmb_path = # path to root directory
# psfriso

aa_translator = {'GLY': 'G',
                 'ALA': 'A',
                 'VAL': 'V',
                 'LEU': 'L',
                 'ILE': 'I',
                 'MET': 'M',
                 'MSE': 'M',
                 'PHE': 'F',
                 'TRP': 'W',
                 'PRO': 'P',
                 'SER': 'S',
                 'THR': 'T',
                 'CYS': 'C',
                 'MEX': 'C',
                 'ABU': 'C',
                 'TYR': 'Y',
                 'ASN': 'N',
                 'GLN': 'Q',
                 'ASP': 'D',
                 'GLU': 'E',
                 'LYS': 'K',
                 'ARG': 'R',
                 'HIS': 'H',
                 'HIE': 'H'}

def existRep(path,trans,i):
    dirpath =path+trans+'-'+str(i)
    return os.path.exists(dirpath) and len(os.listdir(dirpath))>0

def findReplica(path,trans):
    nrep = 0
    for i in xrange(1,100):
        if existRep(path,trans,i):
            nrep = i
            break
    if nrep > 0 :
        return str(nrep)
    else:
        sys.exit('No Replica Found')


def PDBAtomIterator(pdbID,handle):
    #
    def restype(residue):
        return seq1(residue.resname, custom_map=protein_letters_3to1)

    def isX(res):
        return seq1(res.get_resname().upper(), custom_map=protein_letters_3to1 ) == "X"

    undo_handle = UndoHandle( handle )
    struct = PDBParser(QUIET = True).get_structure( pdbID, undo_handle )
    model = struct[0]
    #
    for chn_id, chain in sorted( model.child_dict.items() ):

        residues = [ res for res in chain.get_unpacked_list() if not isX(res) ]

        if not residues: continue

        # Missing residues
        gaps = []
        rnumbers = [ r.id[1] for r in residues ]

        for i,rnum in enumerate( rnumbers[:-1] ):
            if rnumbers[ i+1 ] != rnum +1:
                # its a gap
                gaps.append( (i+1, rnum, rnumbers[ i+1 ] ) )

        if gaps:
            res_out = []
            prev_idx = 0
            # check non-standard aa
            for res in chain.get_unpacked_list():
                if res.get_resname() not in aa_translator:
                    res_out.extend('X')
                    print( res.resname )

            for i, pregap, postgap in gaps:
                if postgap > pregap:
                    gapsize = postgap - pregap - 1
                    res_out.extend( restype(res) for res in residues[prev_idx:i])
                    pred_idx = i
                    res_out.append( 'Z' * gapsize )
                elif postgap == pregap:
                    res_out.extend( restype(res) for res in residues[prev_idx:i])
                    pred_idx = i
                    res_out.append( 'B' )
                else:
                    warnings.warn( "Ignoring out_of_order residues after gap" , BiopythonWarning)
                    res_out.extend( restype(res) for res in residues[ prev_idx : i ])
                    break
            else:
                res_out.extend( restype(res) for res in residues[ prev_idx :  ])

        else:
            # No gaps
            res_out = [ restype(res) for res in residues ]

    #    print( res_out )
        return "".join(res_out)

def get_working_dir(pdb1,pdb2):
    transID = pdb1+'-'+pdb2
    c1 = pdb1[1:3]
    c2 = pdb2[1:3]
    ch1 = pdb1[5]
    ch2 = pdb2[5]
    path = mmb_path+c1+'/'+ch1+'/'+c2+'/'+ch2+'/'
    nrep = findReplica(path,transID)
    return path+transID+'-'+nrep , transID+'-'+nrep

def doIt(fname,pdbfile):
    handle = open(pdbfile,'r')
    f = open( fname , 'w')
    residues = PDBAtomIterator("AAAA" ,handle)
    handle.close()
    #print( residues )
    if 'Z' in residues : f.write("Contains Gaps"+"\n")
    if 'B' in residues : f.write("Contains Insertions"+"\n")
    if 'X' in residues : f.write("Contains non-standard residues"+"\n")
    f.write("Finished PDB Scan"+"\n")
    f.close()


def main(transname):
    pdb1 = transname[0:6]
    pdb2 = transname[7:13]
    trans = pdb1+'-'+pdb2
    wdir, transID = get_working_dir(pdb1,pdb2)
    os.chdir(wdir)
    pdbfile =  mmb_path+pdb1[1:3]+'/'+pdb1+'.pdb'
    doIt('pdb_ini_check.dat',pdbfile)

def overList( fname ):
    f = open(fname,'r')
    for l in f:
        transname = l.rstrip()
        print(transname)
        main(transname)

if __name__ == '__main__':
    overList( sys.argv[1] )
