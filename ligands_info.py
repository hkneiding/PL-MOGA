def get_ligand_names_and_charges(space):

    if space == '1M':
        # --- 1M --- #
        ligand_names = [
            'RUCBEY-subgraph-1', 'WECJIA-subgraph-3', 'KEYRUB-subgraph-1', 'NURKEQ-subgraph-2', 'MEBXUN-subgraph-1',
            'BIFMOV-subgraph-1', 'CUJYEL-subgraph-2', 'EZEXEM-subgraph-1', 'FOMVUB-subgraph-2', 'EFIHEJ-subgraph-3',
            'LETTEL-subgraph-1', 'KAKKIR-subgraph-3', 'BICRIQ-subgraph-3', 'UPEGAZ-subgraph-2', 'CEVJAP-subgraph-2',
            'BABTUT-subgraph-3', 'ZEJJEF-subgraph-3', 'KULGAZ-subgraph-2', 'CIGDAA-subgraph-1', 'HOVMIP-subgraph-3',
            'ULUSIE-subgraph-1', 'IBEKUV-subgraph-1', 'REQSUD-subgraph-2', 'BOSJIF-subgraph-1', 'GUVMEP-subgraph-0',
            'MAZJIJ-subgraph-0', 'OBONEA-subgraph-1', 'CORTOU-subgraph-2', 'LEVGUO-subgraph-2', 'REBWEB-subgraph-2',
            'DOGPAS-subgraph-1', 'IJIMIX-subgraph-1', 'PEJGAN-subgraph-1', 'BIFZEX-subgraph-0', 'IRIXUC-subgraph-3',
            'SAYGOO-subgraph-0', 'UROGIS-subgraph-1', 'MAQKEX-subgraph-1', 'LUQWUQ-subgraph-1', 'QAYDID-subgraph-2',
            'MOYDOV-subgraph-3', 'NIZQUK-subgraph-1', 'SAYHIJ-subgraph-1', 'CIQGOY-subgraph-0', 'VUFZUT-subgraph-1',
            'ZOQFIU-subgraph-0', 'GUQBUQ-subgraph-0', 'LEZYUM-subgraph-2', 'RAJXUX-subgraph-2', 'QEWZOH-subgraph-3'
        ]
        ligand_charges = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
        ]
    elif space == '1B':
        # --- 1B --- #
        ligand_names = [
            'JEZBAR-subgraph-3', 'SASKIG-subgraph-0', 'CEPCAE-subgraph-5', 'NEWYIY-subgraph-5', 'UHIMUU-subgraph-1',
            'SEHHIW-subgraph-1', 'ABESAB-subgraph-0', 'ZUQYEP-subgraph-1', 'USESES-subgraph-0', 'IYAMAX-subgraph-2',
            'WOZRAG-subgraph-1', 'PIWZOM-subgraph-0', 'NEJSUT-subgraph-0', 'PAKLUL-subgraph-2', 'AJUBIS-subgraph-1',
            'PAWGOK-subgraph-1', 'VICTAF-subgraph-1', 'PATPUY-subgraph-2', 'EPOWEM-subgraph-1', 'HOVMIP-subgraph-3',
            'LEFREW-subgraph-0', 'TUNGOZ-subgraph-5', 'ZUBJOV-subgraph-0', 'WODCIB-subgraph-2', 'ERUSEQ-subgraph-2',
            'EGUDAM-subgraph-3', 'RAJVEE-subgraph-2', 'ERURAL-subgraph-3', 'IDOTUT-subgraph-1', 'IGUXAJ-subgraph-2',
            'IWOPOZ-subgraph-1', 'GUBZAG-subgraph-0', 'EFOXII-subgraph-5', 'AGAYOX-subgraph-4', 'YUXWAQ-subgraph-1',
            'OHUTAN-subgraph-3', 'BAHZAN-subgraph-1', 'MULYIC-subgraph-1', 'CBFMOP-subgraph-3', 'XIJWOD-subgraph-1',
            'MOVYAZ-subgraph-3', 'MAYZUH-subgraph-1', 'KEZBEY-subgraph-1', 'NAQWUZ-subgraph-4', 'DEYDIW-subgraph-3',
            'YIMKUA-subgraph-1', 'TIDQAY-subgraph-1', 'SOQHIQ-subgraph-4', 'SICNID-subgraph-3', 'TIJQOS-subgraph-0',
            'YAMTEL-subgraph-4', 'HOVMOV-subgraph-2', 'KATREC-subgraph-3', 'BIFXAP-subgraph-2', 'QEBSAT-subgraph-3',
            'AGUFIT-subgraph-1', 'VEZXAA-subgraph-3', 'RUHPOA-subgraph-1', 'OCEPAQ-subgraph-3', 'KAPZIL-subgraph-2',
            'EZOZIE-subgraph-2', 'VABMET-subgraph-3', 'QUHTUI-subgraph-1', 'TIGXAK-subgraph-4', 'OYOPEY-subgraph-2',
            'HAMTOE-subgraph-3', 'IPELOE-subgraph-3', 'BAXCEK-subgraph-2', 'COPSOQ-subgraph-1', 'ETUVAQ-subgraph-1',
            'KAPDEM-subgraph-1', 'YEXCUA-subgraph-1', 'CECVAK-subgraph-3', 'QETYEU-subgraph-2', 'EZAMOH-subgraph-2',
            'ECUZEI-subgraph-3', 'OBIHOX-subgraph-3', 'RIWYUR-subgraph-1', 'DEBCES-subgraph-2', 'UQUFIW-subgraph-1',
            'IGANOT-subgraph-0', 'LABWUG-subgraph-2', 'EYIWIT-subgraph-2', 'GOMXUB-subgraph-1', 'YESFOR-subgraph-2',
            'CILPET-subgraph-4', 'NOCRED-subgraph-3', 'NORHEJ-subgraph-2', 'WONTOK-subgraph-1', 'KURRIZ-subgraph-4',
            'NEKRAX-subgraph-2', 'YUHNUK-subgraph-2', 'VIDFIZ-subgraph-0', 'NIKGIY-subgraph-1', 'TOPVIE-subgraph-0',
            'DEHDOK-subgraph-1', 'KIJYOS-subgraph-3', 'GOMYAI-subgraph-2', 'QETXET-subgraph-1', 'WUFNAM-subgraph-3',
            'KENTOO-subgraph-1', 'NAXWIT-subgraph-1', 'DUGVUX-subgraph-4', 'CTMOPD-subgraph-3', 'BOFHOU-subgraph-5',
            'QEWGIK-subgraph-1', 'XOCKOS-subgraph-1', 'WAJFOC-subgraph-2', 'GEPGEP-subgraph-3', 'DEJCAX-subgraph-3',
            'QETXOD-subgraph-3', 'CEGZUJ-subgraph-1', 'MEFWAY-subgraph-2', 'BIFZIB-subgraph-3', 'XEXLIV-subgraph-2',
            'ECIBEY-subgraph-0', 'ZOMHAN-subgraph-3', 'JOHDIV-subgraph-1', 'TALXOU-subgraph-1', 'POHKON-subgraph-2',
            'WAJNEC-subgraph-3', 'GOMXOV-subgraph-1', 'TECQOI-subgraph-3', 'PARPEG-subgraph-1', 'NOPXOF-subgraph-3',
            'PELWOT-subgraph-3', 'MAZJIJ-subgraph-0', 'CITFOB-subgraph-0', 'HIRTUX-subgraph-1', 'XILTIV-subgraph-2',
            'YEBFUJ-subgraph-2', 'XECDAN-subgraph-1', 'TEYBOP-subgraph-0', 'PZRMSI-subgraph-1', 'JIJPEX-subgraph-2',
            'CEHNAH-subgraph-2', 'UBOWIR-subgraph-0', 'YOHHIM-subgraph-0', 'SAGFOV-subgraph-1', 'MESWOW-subgraph-3',
            'KIBBAB-subgraph-1', 'ZOPYEK-subgraph-1', 'IFEDED-subgraph-0', 'BAPYEW-subgraph-2', 'COVBEU-subgraph-1',
            'MULRER-subgraph-0', 'VAVROC-subgraph-2', 'KIRVEM-subgraph-3', 'ABODUS-subgraph-0', 'SEBTOI-subgraph-0',
            'BIGYIZ-subgraph-1', 'REPYER-subgraph-1', 'OFIBEK-subgraph-3', 'HUZYUW-subgraph-2', 'KUJHII-subgraph-0',
            'PULFEI-subgraph-2', 'MAKJIU-subgraph-3', 'QIRYOG-subgraph-0', 'FIRCER-subgraph-1', 'HEFTIY-subgraph-1',
            'HUTZOM-subgraph-2', 'EPIMEW-subgraph-1', 'WUWHUT-subgraph-1', 'FIHJAI-subgraph-4', 'BAZTEC-subgraph-3',
            'QAYDOJ-subgraph-2', 'GITWOV-subgraph-0', 'GULYAP-subgraph-3', 'BUBSEW-subgraph-3', 'DAPVUL-subgraph-0',
            'FIJFOW-subgraph-1', 'YUKLIA-subgraph-4', 'NIRKAA-subgraph-1', 'TURSIJ-subgraph-1', 'EMIDEL-subgraph-1',
            'JASVAD-subgraph-0', 'RAFJAJ-subgraph-1', 'HOZHUA-subgraph-0', 'JOBBEK-subgraph-2', 'YOSMUP-subgraph-1',
            'EYICIY-subgraph-2', 'UKOTET-subgraph-1', 'QAVJED-subgraph-3', 'CADJOG-subgraph-2', 'BEWTOO-subgraph-1',
            'ZIDQAG-subgraph-0', 'WEKQAE-subgraph-0', 'ABAJAP-subgraph-0', 'WISJEQ-subgraph-1', 'VAMYUF-subgraph-1',
            'ZAGGIA-subgraph-2', 'ATUROW-subgraph-1', 'NOSNEQ-subgraph-0', 'YACVON-subgraph-2', 'QEHYAD-subgraph-1',
            'XEYQOH-subgraph-1', 'NEWYEV-subgraph-3', 'ASINUL-subgraph-0', 'ABOFAA-subgraph-0', 'QOKDIE-subgraph-4',
            'MORLAH-subgraph-0', 'NACTUH-subgraph-2', 'WECGET-subgraph-1', 'QICLAR-subgraph-1', 'ICEFOL-subgraph-0',
            'BOPSAD-subgraph-2', 'AKARIO-subgraph-1', 'RALMOH-subgraph-0', 'OPEREI-subgraph-1', 'KIFWEB-subgraph-2',
            'RAGDIM-subgraph-1', 'OBOBOY-subgraph-3', 'HONPAC-subgraph-3', 'AREKOY-subgraph-3', 'WAYGAG-subgraph-1',
            'QASDAR-subgraph-1', 'BAPDAY-subgraph-1', 'ITIDUK-subgraph-2', 'BAMXIZ-subgraph-1', 'BAQSER-subgraph-2',
            'KESNUU-subgraph-1', 'HICRAN-subgraph-0', 'YOGWAU-subgraph-0', 'PUMDAC-subgraph-1', 'LAYJIG-subgraph-0',
            'NURLES-subgraph-1', 'MUSPEW-subgraph-1', 'UDOPIM-subgraph-3', 'BAZTOM-subgraph-0', 'UWIWUT-subgraph-0',
            'VEJLEB-subgraph-4', 'FOKKEW-subgraph-2', 'VANTUA-subgraph-1', 'TADREY-subgraph-3', 'YIPZOP-subgraph-4',
            'LOKQAD-subgraph-2', 'XEVDAE-subgraph-0', 'TIXLES-subgraph-3', 'QIVHIO-subgraph-1', 'WOCYET-subgraph-1',
            'EGONEU-subgraph-1', 'JEHKIR-subgraph-5', 'DAJLOS-subgraph-1', 'BAHKOL-subgraph-2', 'FIVSIN-subgraph-0',
            'TOPYEE-subgraph-1', 'SIYVOL-subgraph-0', 'OHEHOZ-subgraph-1', 'TUXGEY-subgraph-0', 'WEMKUV-subgraph-2',
            'KOBVAZ-subgraph-0', 'DACCUG-subgraph-2'
        ]
        ligand_charges = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
        ]
    else:
        raise Exception('Space name not recognized.')
    
    return ligand_names, ligand_charges