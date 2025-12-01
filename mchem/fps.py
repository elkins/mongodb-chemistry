"""
mchem.fps
~~~~~~~~~

Functions for generating fingerprints using RDKit.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from collections import defaultdict
import logging

import pymongo
from rdkit import Chem
from rdkit.Chem import AllChem


log = logging.getLogger(__name__)


def generate(mol_collection, fp_collection, fingerprinter):
    """Generate a fingerprint for all molecules in a collection.

    Example::

        generate(mols, MorganFingerprinter(radius=2))
        generate(mols, MorganFingerprinter(radius=2, length=1024))

    :param mol_collection: MongoDB database collection containing molecules.
    :param fp_collection: MongoDB database collection to store fingerprints.
    :param fingerprinter: fingerprinter instance to generate fingerprint for each molecule.
    """
    log.info(f'Generating {fingerprinter.name} fingerprints for {mol_collection.name} into {fp_collection.name}')
    success, skip = 0, 0
    for molecule in mol_collection.find(no_cursor_timeout=True):
        log.debug(f'Generating {fingerprinter.name} for {molecule["_id"]}')
        bits = fingerprinter.generate(Chem.Mol(molecule['rdmol']))
        fp = {
            '_id': molecule['_id'],
            'bits': bits,
            'count': len(bits)
        }
        try:
            fp_collection.insert_one(fp)
            log.debug(f"Inserted fingerprint for {fp['_id']}")
            success += 1
        except pymongo.errors.DuplicateKeyError:
            log.debug(f"Skipped {fp['_id']}: Fingerprint already exists")
            skip += 1
    log.info(f'{success} successes, {skip} skipped')
    log.info(f'Ensuring index on bits and counts for {fp_collection.name}')
    fp_collection.create_index('bits')
    fp_collection.create_index('count')


def count(fp_collection, count_collection):
    """Build collection containing total counts of all occurrences of each fingerprint bit."""
    counts = defaultdict(int)
    count_collection.drop()
    log.info(f'Counting fingerprint bits in {count_collection.name}')
    for fp in fp_collection.find(no_cursor_timeout=True):
        log.debug(f"Processing {fp['_id']}")
        for bit in fp['bits']:
            counts[bit] += 1
    for k, v in counts.items():
        log.debug(f'Saving count {k}: {v}')
        count_collection.insert_one({'_id': k, 'count': v})


class Fingerprinter(object):
    """Fingerprinter interface."""

    def generate(self, mol):
        """Generate this fingerprint for a molecule."""
        raise NotImplementedError('Fingerprinter subclasses must implement a generate method')

    @property
    def name(self):
        """Unique name for this fingerprint."""
        raise NotImplementedError('Fingerprinter subclasses must implement a name property')


class MorganFingerprinter(Fingerprinter):
    """Class for generating morgan fingerprints."""

    def __init__(self, radius=2, length=None):
        """Initialize with a radius and an optional length.

        :param int radius: The Morgan fingerprint radius (default: 2).
        :param length: The number of bits to optionally fold the fingerprint down to.
        """
        self.radius = radius
        self.length = length

    def generate(self, mol):
        """Generate Morgan fingerprint for a molecule.

        :param mol: The RDKit Mol to generate the fingerprint for.
        """
        if self.length:
            fp = AllChem.GetHashedMorganFingerprint(mol, radius=self.radius, nBits=self.length)
        else:
            fp = AllChem.GetMorganFingerprint(mol, radius=self.radius)
        return sorted(fp.GetNonzeroElements().keys())

    @property
    def name(self):
        """A unique identifier for this fingerprint with the current settings."""
        n = f'm{self.radius}'
        if self.length:
            n = f'{n}l{self.length}'
        return n


def get_fingerprinter(name, radius, length=None):
    fingerprinter = {
        'morgan': MorganFingerprinter(radius=radius, length=length)
        # Add other fingerprinters here in future
    }[name]
    return fingerprinter
