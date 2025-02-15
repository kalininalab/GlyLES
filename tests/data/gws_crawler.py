from pathlib import Path
import requests

import xmltodict


gwb_test_cases= {
    "redEnd--??1D-GalNAc,p(--??1D-Gal,p)--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p)--??1S)--??1S$MONO,Und,-2H,0,redEnd":'Gal(?1-?)GlcNAcOS(?1-?)GalOS(?1-?)GlcNAc(?1-?)[Gal(?1-?)]GalNAc',
    "redEnd--??1D-GalNAc,p(--??1D-Gal,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1S)--??1S$MONO,Und,-2H,0,redEnd":'Gal(?1-?)GlcNAcOS(?1-?)Gal(?1-?)GlcNAcOS(?1-?)Gal(?1-?)GlcNAc(?1-?)[Gal(?1-?)]GalNAc',
    "redEnd--??1D-GalNAc,p(--3b1D-Gal,p(--3a2D-NeuAc,p)--6?1S)--6b1D-GlcNAc,p--4b1D-Gal,p$MONO,Und,-2H,0,redEnd":'Gal(b1-4)GlcNAc(b1-6)[Neu5Ac(a2-3)Gal6S(b1-3)]GalNAc',
    "redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p((--??1D-Gal,p)--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd":'Fuc(?1-?)[Gal(?1-?)]GlcNAcOS(b1-6)[Gal(b1-3)]GalNAc',
    "redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--3b1D-Gal,p--??1D-GlcNAc,p(--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd":'Fuc(?1-?)GlcNAcOS(?1-?)Gal(b1-3)[GlcNAc(b1-6)]GalNAc',
    "redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-Gal,p}--3a2D-NeuAc,p$MONO,Und,-H,0,redEnd":'{Neu5Ac(a2-3)}Gal(b1-3)[Gal(b1-6)]GalNAc',
    "redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p}--??1L-Fuc,p$MONO,Und,-2H,0,redEnd":'{Fuc(a1-?)}Fuc(a1-?)Gal(?1-?)GlcNAc(b1-3)[Fuc(a1-?)Gal(?1-?)GlcNAc(b1-6)]GalNAc',
    "redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--2a1L-Fuc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd":'Fuc(a1-2)Gal(b1-3)[GlcNAc(b1-6)]GalNAc',
    "redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--2a1L-Fuc,p)--6b1D-GlcNAc,p--4b1D-Gal,p(--3a1D-GalNAc,p)--2a1L-Fuc,p$MONO,Und,-2H,0,redEnd":'Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-6)[Fuc(a1-2)Gal(b1-3)]GalNAc',
    "redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--4b1D-Gal,p--4a1D-GlcNAc,p)--6b1D-GlcNAc,p--4b1D-Gal,p--3b1D-GlcNAc,p--4b1D-Gal,p--4a1D-GlcNAc,p$MONO,Und,-2H,0,redEnd":'GlcNAc(a1-4)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-6)[GlcNAc(a1-4)Gal(b1-4)GlcNAc(b1-3)]GalNAc',
    "redEnd--??1D-Glc,p--4b1D-Gal,p(--3a2D-NeuGc,p@270)--4b1D-GalNAc,p$MONO,Und,-H,0,redEnd":'Neu5Gc(a2-3)[GalNAc(b1-4)]Gal(b1-4)Glc',
    "redEnd--??1D-GalNAc,p(--??1D-Gal,p--??1HexNAc,p--??1L-Fuc,p)--??2D-NeuAc,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd":'Fuc(a1-?)HexNAc(?1-?)Gal(?1-?)[Neu5Ac(?2-?)Neu5Ac(?2-?)]GalNAc',
    "redEnd--??1D-GalNAc,p--??2D-KDN,p$MONO,Und,-H,0,redEnd":'Kdn(?2-?)GalNAc',
    "redEnd--?b1D-GlcNAc,p--4b1D-GlcNAc,p--4b1D-Man,p((--3a1D-Man,p--??1D-GlcNAc,p)--4b1D-GlcNAc,p)--6a1D-Man,p--??1D-GlcNAc,p}--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd":'{GlcNAc(?1-?)}GlcNAc(?1-?)Man(a1-3)[GlcNAc(?1-?)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
    "redEnd--?b1D-GlcNAc,p--4b1D-GlcNAc,p--4b1D-Man,p((--3a1D-Man,p--??1D-GlcNAc,p)--6a1D-Man,p--??1D-GlcNAc,p)--4b1D-GlcNAc,p}(--??1D-Gal,p--??2D-NeuAc,p)--??1D-GlcNAc":'{Neu5Ac(?2-?)Gal(?1-?)}{GlcNAc(?1-?)}GlcNAc(?1-?)Man(a1-3)[GlcNAc(?1-?)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
    "redEnd--?b1D-GlcNAc,p--4b1D-GlcNAc,p--4b1D-Man,p((--3a1D-Man,p--2b1D-GlcNAc,p)--4b1D-GlcNAc,p)--6a1D-Man,p--2b1D-GlcNAc,p$MONO,Und,-H,0,redEnd":"GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
}

gwd_files = {
    "https://gitlab.com/glycoinfo/glycoworkbench/-/raw/master/src/main/resources/conf/original_dictionary_files/carbbankraw_dict.gwd?inline=false": "carbbankraw.gwd",
    "https://gitlab.com/glycoinfo/glycoworkbench/-/raw/master/src/main/resources/conf/original_dictionary_files/cfg_dict.gwd?inline=false": "gwb1.gwd",
    "https://gitlab.com/glycoinfo/glycoworkbench/-/raw/master/src/main/resources/conf/original_dictionary_files/glycomedb_dict.gwd?inline=false": "glycomedb.gwd",
    "https://gitlab.com/glycoinfo/glycoworkbench/-/raw/master/src/main/resources/conf/original_dictionary_files/glycosciences_dict.gwd?inline=false": "glycosciences.gwd",
}


def download_file(url, file_name):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(file_name, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download completed.")
    else:
        print(f"Download failed. Status code: {response.status_code}")


def read_gwb_file(file_name):
    data = []
    with open(file_name, "r", encoding="utf-8") as file:
        xml_data = xmltodict.parse(file.read())
    for entry in xml_data["StructureDictionary"]["StructureType"]:
        data.append(entry["@structure"])
    return data


def crawl_gwb():
    root = Path(__file__).parent
    with open(root / "gwb_iupac.tsv", "w") as out:
        for k,v in gwb_test_cases.items():
            out.write(f"{k}\t{v}\n")

    for url, file_name in gwd_files.items():
        download_file(url, root / file_name)


def extract_gwb():
    sequences = set()
    root = Path(__file__).parent
    for file_name in gwd_files.values():
        tmp = set(read_gwb_file(root / file_name))
        sequences = sequences.union(tmp)
    with open(root / "gwb_fuzzy.txt", "w") as out:
        for seq in sequences:
            if len(seq) > 10:
                print(seq, file=out)


if __name__ == "__main__":
    # crawl_gwb()
    extract_gwb()
