"""
mchem.build
~~~~~~~~~~~

Functions for building the database.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

import gzip
import logging
import os
import urllib.request

from bson import Binary
import pymongo
from rdkit import Chem


log = logging.getLogger(__name__)


def download_chembl():
    """Download chembl_19.sdf.gz to data subdirectory."""
    url = 'ftp://ftp.ebi.ac.uk//pub/databases/chembl/ChEMBLdb/releases/chembl_19/chembl_19.sdf.gz'
    dest = os.path.realpath(os.path.join(os.path.dirname(__file__), '../data/chembl_19.sdf.gz'))
    log.info(f'Downloading {url}')
    if not os.path.isfile(dest):
        r = urllib.request.urlopen(url)
        with open(dest, 'wb') as f:
            f.write(r.read())
    else:
        log.info(f'Already downloaded: {dest}')


def load_sdfs(paths, collection, idfield=None):
    """Insert molecules from list of SDF file paths into database collection.

    :param string paths: Paths to SDF files.
    :param collection: MongoDB collection.
    :param string idfield: SDF property field to use for molecule ID.
    """
    log.info(f'Inserting molecules into MongoDB {collection.name}')
    for path in paths:
        log.info(f'Loading {path}')
        sdf = gzip.open(path, 'rt') if path.endswith('.gz') else open(path)
        load_sdf(sdf, collection, idfield)


def load_sdf(sdf, collection, idfield):
    """Insert molecules from SDF file into database collection.

    :param file sdf: SDF file object.
    :param collection: MongoDB collection.
    :param string idfield: SDF property field to use for molecule ID.
    """

    success, fail, skip = 0, 0, 0
    for rdmol in Chem.ForwardSDMolSupplier(sdf):
        if rdmol is None:
            log.debug('Failed to read molecule')
            fail += 1
            continue
        try:
            smiles = Chem.MolToSmiles(rdmol, isomericSmiles=True)
        except (ValueError, RuntimeError):
            log.debug('Failed to generate SMILES')
            fail += 1
            continue
        mol = {
            'smiles': smiles,
            'rdmol': Binary(rdmol.ToBinary()),
        }
        if idfield:
            mol['_id'] = rdmol.GetProp(idfield)
        try:
            collection.insert_one(mol)
            log.debug(f"Inserted {mol['_id']}")
            success += 1
        except pymongo.errors.DuplicateKeyError:
            log.debug(f"Skipped {mol['_id']}: Already exists")
            skip += 1
    log.info(f'{success} successes, {fail} failures, {skip} skipped')


