#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
Guide Mission 1 - Bioinformatique et TP
=============================================================================

Ce script reproduit le notebook Google Colab "Guided Mission 1" pour VS Code.
Il effectue l'analyse de données de séquençage (RNA-seq, CLIP-seq, RPF-seq)
pour étudier l'effet de Lin28a sur la traduction.

Auteur: Adaptation du notebook original de 장혜식 (Hyeshik Chang)
Cours: 생물정보학 및 실습 1 (Bioinformatique et TP 1)
Institution: Seoul National University

=============================================================================
DESCRIPTION DES DONNÉES
=============================================================================

Le projet utilise plusieurs fichiers BAM (Binary Alignment Map):
- CLIP-35L33G.bam: Données CLIP-seq (crosslinking immunoprecipitation)
- RNA-control.bam: RNA-seq de contrôle
- RNA-siLin28a.bam: RNA-seq après knockdown de Lin28a
- RNA-siLuc.bam: RNA-seq avec knockdown Luciférase (contrôle)
- RPF-siLin28a.bam: Ribosome footprint (RPF) après knockdown Lin28a
- RPF-siLuc.bam: Ribosome footprint contrôle

L'objectif est de:
1. Compter les reads par gène avec featureCounts
2. Calculer l'enrichissement CLIP (CLIP/RNA-control)
3. Calculer le changement de densité ribosomique (rden_change)
4. Créer des scatter plots similaires aux figures de l'article

=============================================================================
"""

import os
import sys
import subprocess
import ssl
import urllib.request
import hashlib

# ============================================================================
# CONFIGURATION - MODIFIEZ CES VARIABLES SELON VOTRE ENVIRONNEMENT
# ============================================================================

# Répertoire de travail (où les données seront stockées)
WORK_DIR = os.path.expanduser("~/binfo_project")
DATA_PACK_DIR = os.path.join(WORK_DIR, "binfo1-datapack1")
WORK_SUBDIR = os.path.join(WORK_DIR, "binfo1-work")

# URLs des données
DATA_URL = "https://hyeshik.qbio.io/binfo/binfo1-datapack1.tar"
GENCODE_URL = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.annotation.gtf.gz"
LOCALIZATION_URL = "https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt"

# MD5 checksums pour vérification d'intégrité
EXPECTED_MD5 = {
    "CLIP-35L33G.bam": "140aaf30bcb9276cc716f8699f04ddd6",
    "CLIP-35L33G.bam.bai": "f1b3336ed7e2f97d562dcc71641251bd",
    "RNA-control.bam": "328883a73d507eafbf5b60bd6b906201",
    "RNA-control.bam.bai": "02073818e2f398a73c3b76e5169de1ca",
    "RNA-siLin28a.bam": "b09550d09d6c2a4ce27f0226f426fdb1",
    "RNA-siLin28a.bam.bai": "fef112c727244060ea62d3f2564a07f6",
    "RNA-siLuc.bam": "28bbd0c47d725669340c784f1b772c01",
    "RNA-siLuc.bam.bai": "43590fdc4d81905c0432e0d1cb8cfd5b",
    "RPF-siLin28a.bam": "5c08a9297307bc83259e658c4474f0cc",
    "RPF-siLin28a.bam.bai": "a1bb3e29be412dfd7fd8d16b1b1acc4c",
    "RPF-siLuc.bam": "f2eebf50943024d0116c9cd3e744c707",
    "RPF-siLuc.bam.bai": "dc24f69e8f571fc8be30f28ce5b84fcd",
}

# ============================================================================
# FONCTIONS UTILITAIRES
# ============================================================================

def run_command(cmd, description="", shell=True):
    """
    Exécute une commande shell et affiche le résultat.
    
    Args:
        cmd: Commande à exécuter (string ou list)
        description: Description de l'opération
        shell: Si True, exécuter comme commande shell
    
    Returns:
        tuple: (success, output)
    """
    print(f"\n{'='*60}")
    if description:
        print(f"  {description}")
    print(f"{'='*60}")
    print(f"Commande: {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(
            cmd if isinstance(cmd, str) else cmd,
            shell=shell,
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        print(f"ERREUR: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False, e.stderr
    except Exception as e:
        print(f"ERREUR INATTENDUE: {e}")
        return False, str(e)


def download_file(url, output_path, description="Téléchargement"):
    """
    Télécharge un fichier depuis une URL.
    
    Args:
        url: URL du fichier à télécharger
        output_path: Chemin de destination
        description: Description pour l'affichage
    
    Returns:
        bool: True si succès
    """
    print(f"\n{description}...")
    print(f"URL: {url}")
    print(f"Destination: {output_path}")
    
    # Créer le contexte SSL qui ignore les erreurs de certificat
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"✓ Téléchargement terminé: {output_path}")
        return True
    except Exception as e:
        print(f"✗ Erreur lors du téléchargement: {e}")
        return False


def verify_md5(file_path, expected_md5):
    """
    Vérifie le checksum MD5 d'un fichier.
    
    Args:
        file_path: Chemin du fichier
        expected_md5: Somme MD5 attendue
    
    Returns:
        bool: True si le checksum correspond
    """
    print(f"Vérification MD5 de: {os.path.basename(file_path)}")
    md5_hash = hashlib.md5()
    
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    
    actual_md5 = md5_hash.hexdigest()
    match = actual_md5 == expected_md5
    
    if match:
        print(f"  ✓ MD5 vérifié: {actual_md5}")
    else:
        print(f"  ✗ ERREUR MD5!")
        print(f"    Attendu: {expected_md5}")
        print(f"    Obtenu: {actual_md5}")
    
    return match


def create_directory(path):
    """
    Crée un répertoire s'il n'existe pas.
    """
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Répertoire créé: {path}")
    else:
        print(f"Répertoire existant: {path}")


# ============================================================================
# ÉTAPE 1: INSTALLATION DES PACKAGES PYTHON
# ============================================================================

def step1_install_python_packages():
    """
    Installe les packages Python nécessaires avec pip.
    
    Packages installés:
    - pandas: Manipulation de données tabulaires
    - matplotlib: Visualisation de graphiques
    - numpy: Calculs numériques
    - scipy: Calculs scientifiques supplémentaires
    """
    print("\n" + "="*70)
    print("  ÉTAPE 1: Installation des packages Python")
    print("="*70)
    
    packages = ["pandas", "matplotlib", "numpy", "scipy"]
    
    for package in packages:
        print(f"\nInstallation de {package}...")
        success, _ = run_command(
            f"pip install {package}",
            f"Installation de {package}"
        )
        if success:
            print(f"✓ {package} installé avec succès")
        else:
            print(f"⚠ {package} pourrait déjà être installé ou erreur")
    
    print("\n✓ Étape 1 terminée: Packages Python installés")


# ============================================================================
# ÉTAPE 2: INSTALLATION DE BIOCONDA ET SUBREAD
# ============================================================================

def step2_install_bioconda():
    """
    Installe Bioconda et le package subread (contenant featureCounts).
    
    Bioconda est un gestionnaire de packages pour la bioinformatique.
    featureCounts est un outil pour compter les reads alignés sur les gènes.
    """
    print("\n" + "="*70)
    print("  ÉTAPE 2: Installation de Bioconda et subread")
    print("="*70)
    
    # Vérifier si conda est déjà installé
    conda_check, _ = run_command("which conda", "Vérification de conda")
    
    if conda_check:
        print("Conda déjà installé!")
    else:
        print("Installation de Miniconda...")
        # Installer Miniconda silencieusement
        run_command(
            "curl -sL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o ~/miniconda.sh",
            "Téléchargement de Miniconda"
        )
        run_command(
            "bash ~/miniconda.sh -b -p ~/miniconda3",
            "Installation de Miniconda"
        )
        
        # Ajouter au PATH
        os.environ["PATH"] = os.path.expanduser("~/miniconda3/bin") + ":" + os.environ["PATH"]
    
    # Installer subread via conda
    print("\nInstallation de subread (featureCounts)...")
    success, _ = run_command(
        "conda install -y -c bioconda subread",
        "Installation de subread via Bioconda"
    )
    
    if success:
        print("✓ subread installé avec succès")
    else:
        print("⚠ Erreur possible lors de l'installation de subread")
    
    print("\n✓ Étape 2 terminée: Bioconda et subread installés")


# ============================================================================
# ÉTAPE 3: TÉLÉCHARGEMENT DES DONNÉES
# ============================================================================

def step3_download_data():
    """
    Télécharge les données de séquençage depuis le serveur.
    
    Cette étape télécharge:
    - 6 fichiers BAM (données d'alignement)
    - 6 fichiers BAI (index des fichiers BAM)
    - 1 fichier GTF (annotations GENCODE)
    
    Taille totale: environ 2 GB
    """
    print("\n" + "="*70)
    print("  ÉTAPE 3: Téléchargement des données de séquençage")
    print("="*70)
    
    # Créer les répertoires
    create_directory(WORK_DIR)
    create_directory(DATA_PACK_DIR)
    create_directory(WORK_SUBDIR)
    
    # Télécharger l'archive des données
    print("\n--- Téléchargement des fichiers BAM et BAI ---")
    archive_path = os.path.join(WORK_DIR, "binfo1-datapack1.tar")
    
    success = download_file(DATA_URL, archive_path, "Archive des données")
    
    if success:
        # Extraire l'archive
        print("\nExtraction de l'archive...")
        run_command(
            f"tar -xf {archive_path} -C {WORK_DIR}",
            "Extraction des fichiers"
        )
        print(f"✓ Données extraites dans: {DATA_PACK_DIR}")
    
    # Télécharger GENCODE GTF
    print("\n--- Téléchargement des annotations GENCODE ---")
    gencode_gz = os.path.join(DATA_PACK_DIR, "gencode.gtf.gz")
    success = download_file(GENCODE_URL, gencode_gz, "Annotations GENCODE")
    
    if success:
        # Décompresser
        run_command(
            f"gunzip -f {gencode_gz}",
            "Décompression de GENCODE"
        )
        print("✓ GENCODE décompressé")
    
    # Copier les données vers le répertoire de travail
    print("\n--- Copie des données vers le répertoire de travail ---")
    run_command(
        f"cp -r {DATA_PACK_DIR}/* {WORK_SUBDIR}/",
        "Copie des fichiers"
    )
    
    print("\n✓ Étape 3 terminée: Données téléchargées et préparées")


# ============================================================================
# ÉTAPE 4: VÉRIFICATION D'INTÉGRITÉ
# ============================================================================

def step4_verify_data():
    """
    Vérifie l'intégrité des fichiers téléchargés avec MD5.
    
    Chaque fichier doit correspondre au checksum attendu pour garantir
    qu'il n'y a pas eu d'erreur de téléchargement ou de corruption.
    """
    print("\n" + "="*70)
    print("  ÉTAPE 4: Vérification d'intégrité (MD5)")
    print("="*70)
    
    all_valid = True
    
    for filename, expected_md5 in EXPECTED_MD5.items():
        file_path = os.path.join(WORK_SUBDIR, filename)
        
        if os.path.exists(file_path):
            if not verify_md5(file_path, expected_md5):
                all_valid = False
        else:
            print(f"✗ Fichier manquant: {filename}")
            all_valid = False
    
    if all_valid:
        print("\n✓ Tous les fichiers sont intègres!")
    else:
        print("\n⚠ Certains fichiers ont des problèmes d'intégrité")
    
    print("\n✓ Étape 4 terminée: Vérification d'intégrité effectuée")


# ============================================================================
# ÉTAPE 5: EXÉCUTION DE FEATURECOUNTS
# ============================================================================

def step5_run_featurecounts():
    """
    Exécute featureCounts pour compter les reads par gène.
    
    featureCounts est un outil qui compte le nombre de reads alignés
    sur chaque feature (gène, exon, transcript) défini dans le fichier GTF.
    
    Input:
    - gencode.gtf: Fichier d'annotation GENCODE (annotations des gènes)
    - *.bam: Fichiers d'alignement (6 fichiers)
    
    Output:
    - read-counts.txt: Tableau avec les counts par gène
    """
    print("\n" + "="*70)
    print("  ÉTAPE 5: Exécution de featureCounts")
    print("="*70)
    
    # Se déplacer dans le répertoire de travail
    os.chdir(WORK_SUBDIR)
    print(f"Répertoire de travail: {os.getcwd()}")
    
    # Vérifier les fichiers disponibles
    print("\nFichiers disponibles:")
    run_command("ls -lh *.bam *.gtf", "Liste des fichiers")
    
    # Exécuter featureCounts
    # -a: fichier d'annotation
    # -o: fichier de sortie
    # *.bam: tous les fichiers BAM
    print("\nExécution de featureCounts...")
    success, output = run_command(
        "featureCounts -a gencode.gtf -o read-counts.txt *.bam",
        "Comptage des reads par gène"
    )
    
    if success:
        print("✓ featureCounts exécuté avec succès")
    else:
        print("⚠ Erreur lors de l'exécution de featureCounts")
    
    # Vérifier le fichier de sortie
    if os.path.exists("read-counts.txt"):
        print(f"\nFichier de sortie créé: read-counts.txt")
        run_command("head -5 read-counts.txt", "Aperçu du résultat")
    
    print("\n✓ Étape 5 terminée: Counts générés")


# ============================================================================
# ÉTAPE 6: ANALYSE DES DONNÉES
# ============================================================================

def step6_analyze_data():
    """
    Analyse les données et crée des visualisations.
    
    Cette étape:
    1. Charge les counts dans pandas
    2. Calcule l'enrichissement CLIP (CLIP/RNA-control)
    3. Calcule le changement de densité ribosomique
    4. Crée un scatter plot basique
    """
    print("\n" + "="*70)
    print("  ÉTAPE 6: Analyse des données et visualisation")
    print("="*70)
    
    # Import des bibliothèques
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Charger les données de counts
    print("\nChargement des données de counts...")
    cnts = pd.read_csv(
        os.path.join(WORK_SUBDIR, 'read-counts.txt'),
        sep='\t',
        comment='#',
        index_col=0
    )
    
    print(f"Nombre de gènes: {len(cnts)}")
    print(f"Colonnes: {list(cnts.columns)}")
    print("\nAperçu des données:")
    print(cnts.head())
    
    # Calculer les métriques
    print("\n--- Calcul des métriques ---")
    
    # Enrichissement CLIP: ratio CLIP/RNA-control
    # Indique quels gènes sont enrichis dans les données CLIP
    cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
    
    # Changement de densité ribosomique:
    # (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc)
    # Compare l'efficacité de traduction entre knockdown Lin28a et contrôle
    cnts['rden_change'] = (
        cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']
    ) / (
        cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam']
    )
    
    print("\nMétriques calculées:")
    print(cnts[['clip_enrichment', 'rden_change']].head(10))
    
    # Créer le scatter plot basique
    print("\n--- Création du scatter plot ---")
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    
    # Scatter plot avec échelle log2
    ax.scatter(
        np.log2(cnts['clip_enrichment']),
        np.log2(cnts['rden_change']),
        alpha=0.5,
        s=5
    )
    
    ax.set_xlabel('log2(CLIP Enrichment)')
    ax.set_ylabel('log2(Ribosome Density Change)')
    ax.set_title('CLIP Enrichment vs Ribosome Density Change')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    # Sauvegarder le graphique
    output_path = os.path.join(WORK_SUBDIR, 'scatter_basic.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Graphique sauvegardé: {output_path}")
    plt.close()
    
    print("\n✓ Étape 6 terminée: Analyse et visualisation basique effectuées")


# ============================================================================
# ÉTAPE 7: GRAPHIQUE AVEC LOCALISATION PROTÉIQUE
# ============================================================================

def step7_localization_plot():
    """
    Crée un scatter plot coloré selon la localisation cellulaire.
    
    Cette étape reproduit les figures Figure 5B et S6A de l'article,
    où les points sont colorés selon la localisation subcellulaire
    des protéines (cytoplasme, noyau, mitochondrie, etc.)
    """
    print("\n" + "="*70)
    print("  ÉTAPE 7: Graphique avec localisation protéique")
    print("="*70)
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Charger les données de counts
    cnts = pd.read_csv(
        os.path.join(WORK_SUBDIR, 'read-counts.txt'),
        sep='\t',
        comment='#',
        index_col=0
    )
    
    # Recalculer les métriques
    cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
    cnts['rden_change'] = (
        cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']
    ) / (
        cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam']
    )
    
    # Télécharger les données de localisation
    print("\nTéléchargement des données de localisation...")
    
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    
    local_path = os.path.join(WORK_SUBDIR, 'mouselocalization.txt')
    urllib.request.urlretrieve(LOCALIZATION_URL, local_path)
    
    mouselocal = pd.read_csv(local_path, sep='\t')
    print(f"Données de localisation chargées: {len(mouselocal)} entrées")
    print(mouselocal.head())
    
    # Préparer les noms de gènes (enlever la version)
    cnts['gene_id'] = cnts.index.str.split('.').str[0]
    mouselocal['gene_id'] = mouselocal['Gene ID'].str.split('.').str[0]
    
    # Fusionner les données
    merged = cnts.merge(
        mouselocal[['gene_id', 'Localizations']],
        on='gene_id',
        how='left'
    )
    
    print(f"\nDonnées fusionnées: {len(merged)} gènes")
    print(f"Gènes avec localisation: {merged['Localizations'].notna().sum()}")
    
    # Créer le graphique avec couleurs par localisation
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Définir les couleurs pour les localisations principales
    localization_colors = {
        'Cytoplasm': '#1f77b4',      # Bleu
        'Nucleus': '#ff7f0e',        # Orange
        'Mitochondrion': '#2ca02c',  # Vert
        'Endoplasmic reticulum': '#d62728',  # Rouge
        'Secreted': '#9467bd',       # Violet
    }
    
    # Tracer les points sans localisation (gris)
    no_local = merged['Localizations'].isna()
    ax.scatter(
        np.log2(merged.loc[no_local, 'clip_enrichment']),
        np.log2(merged.loc[no_local, 'rden_change']),
        c='lightgray', alpha=0.3, s=5, label='Unknown'
    )
    
    # Tracer les points avec localisation
    has_local = ~no_local
    for loc, color in localization_colors.items():
        mask = merged['Localizations'].str.contains(loc, na=False)
        if mask.any():
            ax.scatter(
                np.log2(merged.loc[mask, 'clip_enrichment']),
                np.log2(merged.loc[mask, 'rden_change']),
                c=color, alpha=0.6, s=15, label=loc
            )
    
    ax.set_xlabel('log2(CLIP Enrichment)', fontsize=12)
    ax.set_ylabel('log2(Ribosome Density Change)', fontsize=12)
    ax.set_title('CLIP Enrichment vs Ribosome Density Change\n(Colored by Protein Localization)', fontsize=14)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize=10)
    
    # Sauvegarder
    output_path = os.path.join(WORK_SUBDIR, 'scatter_with_localization.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nGraphique sauvegardé: {output_path}")
    plt.close()
    
    print("\n✓ Étape 7 terminée: Graphique avec localisation créé")


# ============================================================================
# FONCTION PRINCIPALE
# ============================================================================

def main():
    """
    Fonction principale qui exécute toutes les étapes.
    
    Vous pouvez commenter/décommenter les étapes selon vos besoins.
    """
    print("\n" + "="*70)
    print("  GUIDE MISSION 1 - BIOINFORMATIQUE")
    print("  Adaptation pour VS Code")
    print("="*70)
    
    print(f"\nRépertoire de travail: {WORK_DIR}")
    print(f"Répertoire des données: {DATA_PACK_DIR}")
    print(f"Répertoire de travail: {WORK_SUBDIR}")
    
    # Décommenter les étapes que vous voulez exécuter:
    
    # Étape 1: Installer les packages Python
    step1_install_python_packages()
    
    # Étape 2: Installer Bioconda et subread
    step2_install_bioconda()
    
    # Étape 3: Télécharger les données
    step3_download_data()
    
    # Étape 4: Vérifier l'intégrité
    step4_verify_data()
    
    # Étape 5: Exécuter featureCounts
    step5_run_featurecounts()
    
    # Étape 6: Analyse basique
    step6_analyze_data()
    
    # Étape 7: Graphique avec localisation
    step7_localization_plot()
    
    print("\n" + "="*70)
    print("  TOUTES LES ÉTAPES TERMINÉES!")
    print("="*70)
    print(f"\nRésultats sauvegardés dans: {WORK_SUBDIR}")
    print("\nFichiers générés:")
    print("  - read-counts.txt: Tableau des counts par gène")
    print("  - scatter_basic.png: Graphique basique")
    print("  - scatter_with_localization.png: Graphique avec localisation")


if __name__ == "__main__":
    main()