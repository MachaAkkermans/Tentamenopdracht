import matplotlib.pyplot as plt
import re
import tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


class Gff_bestand():

    def __init__(self, bestandsnaam):
        self.__bestandsnaam = bestandsnaam
        self.bestand_in_lijst()
        self.lijsten_aanmaken()

    def bestand_in_lijst(self):
        """
        Maakt lijst self.__inhoud_gff aan, waar elke regel van het bestand
        wordt toegevoegt.
        """
        try:
            bestand = open(self.__bestandsnaam)
            self.__inhoud_gff = []
            for regel in bestand:
                regel = regel.replace("\n", "")  # Haalt de enter eruit
                lijst = regel.split("\t")  # splits de lijst op elke tab
                self.__inhoud_gff.append(lijst)
            bestand.close()
        except FileNotFoundError:
            print("bestand niet gevonden!")
            self.__bestandsnaam = input("voer bestand opnieuw in")
            self.__init__(self.__bestandsnaam)
            self.getlijsten()

    def lijsten_aanmaken(self):
        """
        Maakt 5 lijsten aan, voor elke soort type die voorkomt.
        """
        self.__cds = []
        self.__gene = []
        self.__trna = []
        self.__exon = []
        self.__anders = []
        # stopt elke regel van het bestand in de juiste lijst
        for lijst in self.__inhoud_gff:
            if lijst[2] == "gene":
                self.__gene.append(lijst)
            elif lijst[2] == "CDS":
                self.__cds.append(lijst)
            elif lijst[2] == "tRNA":
                self.__trna.append(lijst)
            elif lijst[2] == "exon":
                self.__exon.append(lijst)
            else:
                self.__anders.append(lijst)

    def getlijsten(self):
        """

        :return:
        self.__cds: lijst met regels die informatie bevatten over het cds
        self.__gene: lijst met regels die informatie bevatten over het gene
        self.__trna: lijst met regels die informatie bevatten over en tRNA
        self.__exon: lijst met regels die informatie bevatten over het exon
        self.__anders: lijst met regels die informatie bevatten over de rest
        zoals, region en pseudogene
        self.__inhoud_gff: 2d lijst waarin elke regel van het bestand als
        lijst in de 2d lijst is toegevoegt:
        [[NC_007795.1, Refseq, gene, ...],[NC_007795.1, Refseq, CDS,...] N]
        """
        return self.__cds, self.__gene, self.__trna, self.__exon, \
               self.__anders, self.__inhoud_gff


class Gffb_bestand():
    def __init__(self, bestandsnaam):
        self.__bestandsnaam = bestandsnaam
        self.bestand_opdelen()
        self.informatie_halen_uit_bestand2()
        self.dictonaries_aanmaken()
        # self.sequentie_uit_bestand3()

    def bestand_opdelen(self):
        """
        splits het bestand op in 3 bestanden
        # bestand 1 bevat informatie van LOCUS t/m FEATURES
        # bestand 2 bevat informatie over de CDS, gene enz.
        # bestand 3 bevat de sequentie
        :return:
        """
        self.__bestand1_tekst = open("bestand1", "w")
        self.__bestand2_informatie = open("gene+cds.txt", "w")
        self.__bestand3_sequentie = open("sequentie.txt", "w")
        try:
            bestand = open(self.__bestandsnaam, "r")

            # houdt bij bij welke van de 3 bestanden het de regel moet toevoegen

            bestand1 = True
            bestand2 = False
            bestand3 = False
            for regel in bestand:
                if regel.startswith("FEATURES"):  # vanaf hier begint deel 2
                    # dus wordt bestand1 op false gezet en bestand 2 op True
                    bestand1 = False
                    bestand2 = True
                if regel.startswith("CONTIG"):  # vanaf hier begint deel 3
                    # dus wordt bestand 2 op false gezet en bestand3 op true
                    bestand2 = False
                    bestand3 = True
                # de code hieronder zorgt dat de regel bij het goede bestand
                # wordt toegevoegt
                if bestand1 == True and bestand2 == False:
                    self.__bestand1_tekst.write(regel)
                elif bestand1 == False and bestand2 == True \
                        and bestand3 == False:
                    self.__bestand2_informatie.write(regel)
                elif bestand2 == False and bestand3 == True:
                    self.__bestand3_sequentie.write(regel)

            self.__bestand1_tekst.close()
            self.__bestand2_informatie.close()
            self.__bestand3_sequentie.close()
        except FileNotFoundError:
            print("bestand niet gevonden!")
            self.__bestandsnaam = input("voer bestand opnieuw in")
            self.__init__(self.__bestandsnaam)
            self.getdiconaries()
            self.getsequentie()

    def spaties_uit_regel(self):
        """
        voorbeeld van recursie om spaties in regel eruit te halen
        """
        if self.__regel.startswith(" "):
            self.__regel = self.__regel[1:]
            self.spaties_uit_regel()

    def informatie_halen_uit_bestand2(self):
        """

        :return: lijst_met_gene_cds - 2d lijst met inhoudt van bestand gene+cds.txt
        - [["gene 517..1878", /gene="dnaA", ...] [CDS 517..1878, /gene="dnaA",...],
          ["gene ..., n]]
        """
        bestand = open(self.__bestand2_informatie.name)
        # bestand = open("test.txt")
        self.__lijst_met_gegevens = []
        per_gegeven = []
        per_kopje_in_bestand = ""
        for self.__regel in bestand:
            if self.__regel.startswith(
                    " "):  # Haalt alle spaties weg dat bij het
                self.spaties_uit_regel()  # begin van de regel staat
            self.__regel = self.__regel.replace("\n", "")
            # zorgt ervoor dat elke keer wanneer "gene" of "CDS"staat er een nieuw
            # lijst wordt aangemaakt en de oude wordt toegevoegt aan de 2d lijst
            if self.__regel.startswith("gene") or self.__regel.startswith(
                    "CDS") or self.__regel.startswith("rRNA   "):
                per_gegeven.append(per_kopje_in_bestand)
                self.__lijst_met_gegevens.append(per_gegeven)
                per_gegeven = []
                per_gegeven.append(self.__regel)
                per_kopje_in_bestand = ""
            else:
                # Zorgt ervoor dat kopjes die bestaan uit meerdere regels
                # samen als een regel worden toegevoegt in een lijst
                if self.__regel.startswith("/"):
                    if per_kopje_in_bestand != "":
                        per_gegeven.append(per_kopje_in_bestand)
                        per_kopje_in_bestand = ""
                    per_kopje_in_bestand += self.__regel
                else:
                    per_kopje_in_bestand += self.__regel
        per_gegeven.append(per_kopje_in_bestand)
        self.__lijst_met_gegevens.append(per_gegeven)
        del self.__lijst_met_gegevens[0]  # verwijdert de eerste legen
        # inhoudt

    def dictonaries_aanmaken(self):
        """
        Maakt 3 dictonaries aan, elk voor elk type dat voorkomt in het bestand
        """
        self.__gene = {}
        self.__cds = {}
        self.__rrna = {}
        # houdt ID_nummers bij om te gebruiken als key in de dictonary
        id_nummer_gene = 1
        id_nummer_cds = 1
        id_nummer_rrna = 1
        for lijst in self.__lijst_met_gegevens:
            # Kijkt in welke dictonary het hoort en voegt en daar aan toe
            if lijst[0].startswith("gene  "):
                self.__gene[id_nummer_gene] = lijst
                id_nummer_gene += 1
            elif lijst[0].startswith("CDS  "):
                self.__cds[id_nummer_cds] = lijst
                id_nummer_cds += 1
            else:
                self.__rrna[id_nummer_rrna] = lijst
                id_nummer_rrna += 1

    def getsequentie(self):
        """

        :return: self.__bestand3__sequentie - naam van bestand 3
        """
        return self.__bestand3_sequentie

    def getdiconaries(self):
        """

        :return:
        self.__gene - diconary met gegevens over genen
        { 1 : [gene  2156..3289,/locus_tag="SAOUHSC_00002",..]
          2: n
        self.__cds - diconary met gegevens over CDS
        { 1: [CDS   2156..3289,/locus_tag="SAOUHSC_00002",..]
          2: n
        self.__rrna - diconary met gegevens over rrna
        { 1: [rRNA    448819..450374,/locus_tag="SAOUHSC_R0001,..]
          2: n
        """
        return self.__gene, self.__cds, self.__rrna


class Gui_maken():
    def __init__(self, patronen, gff_totale_inhoud, gffb_cds, sequentie,
                 gff_gene):
        self.__patronen = patronen
        self.__gff_totale_inhoud = gff_totale_inhoud
        self.__gffb_cds = gffb_cds
        self.__sequentie = sequentie
        self.__gff_gene = gff_gene
        self.setgegevenstypes()
        self.setgegevenshypothetical_proteins()
        self.settabel()
        self.setgegevens_strands()
        self.setgegevens_gc_percentage()
        self.setLengte()

        self.gui_startscherm()

    def gui_startscherm(self):
        """
        Bouwt de GUI en knopjes voor de verschillende commands
        """
        startscherm = Tk.Tk()
        startscherm.title("informatie over Staphylococcus aureus")
        self.__knop = Tk.Button(heigh=2, master=startscherm,
                                text="klik hier voor tabel met de "
                                     "gevonden patronen",
                                command=self.tabel_patronen)
        self.__knop2 = Tk.Button(heigh=2, master=startscherm,
                                 text="klik hier voor grafiek met "
                                      "verschillende gevonden type ",
                                 command=self.grafiek_type)
        self.__knop3 = Tk.Button(heigh=2, master=startscherm, text=
                                "klik hier voor grafiek met het aantal "
                                    "hypotheticale proteins ",
                                 command=
                                 self.grafiek_hypothetical_proteins)
        self.__knop4 = Tk.Button(heigh=2, master=startscherm,
                                 text="klik hier voor grafiek met het "
                                      "verschil +/- strands ",
                                 command=self.maak_grafiek_strands)
        self.__knop5 = Tk.Button(heigh=2, master=startscherm, text=
                                 "klik hier voor grafiek met GC/AT %",
                                 command=self.maak_grafiek_percentage)
        self.__knop6 = Tk.Button(heigh=2, master=startscherm, text=
                         "klik hier verschillende lengtes van genen",
                                 command=self.maakGrafiek_lengte)
        quit_button = Tk.Button(heigh=2, width=10, text="quit",
                                command=startscherm.destroy, )

        self.__knop.pack()
        self.__knop2.pack()
        self.__knop3.pack()
        self.__knop4.pack()
        self.__knop5.pack()
        self.__knop6.pack()
        quit_button.pack()
        Tk.mainloop()

    def setgegevenstypes(self):
        """
        Berekent de gegevens voor de grafiek verschillende voorkomende type
        """
        self.__labels = {}  # dictonary met alle unieke labels en hoevaak ze
        # voor komen
        for lijst in self.__gff_totale_inhoud:
            try:
                if lijst[2] not in self.__labels:
                    self.__labels[lijst[2]] = 1
                else:
                    self.__labels[lijst[2]] += 1
            except IndexError:
                print("List out of range!")
        self.__x = []  # waarde wat op de x-as moet komen te staan
        self.__y = []  # waarde wat op de y-as moet komen te staan
        for key in self.__labels.keys():
            self.__x.append(key)
        for value in self.__labels.values():
            self.__y.append(value)

    def settabel(self):
        """
        Maakt een lijst aan met de gegevens wat in de tabel moet komen
        te staan
        """
        self.__tabel = []
        # Gegevens voor de kopjes
        self.__tabel.append(
            ["nummer", "eiwit", "protein_id", "product", "+/- strand"])
        # Gegevens voor de inhoud
        for patroon in self.__patronen:
            self.__tabel.append(patroon)

    def setgegevens_gc_percentage(self):
        """
        Berekent de gegevens voor de pie chart
        """
        try:
            self.__gc = round(
                (self.__sequentie.count("c") + self.__sequentie.count
                ("g")) / len(self.__sequentie) * 100, 2)
            # Berekent het GC% en rond af op 2
        except ZeroDivisionError:
            self.__gc = 0
        try:
            self.__at = round(
                (self.__sequentie.count("a") + self.__sequentie.count
                ("t")) / len(self.__sequentie) * 100, 2)
            # Berekent het AT% en rond af op 2
        except ZeroDivisionError:
            self.__at = 0
        self.__percentage = [self.__gc, self.__at]

    def setgegevenshypothetical_proteins(self):
        """
        Berekent de aantal hypothetical proteins
        """
        aantal_wel = 0
        aantal_niet = 0
        for k, v in self.__gffb_cds.items():
            for item in v:
                if item == "/note=\"conserved hypothetical protein\"":
                    aantal_wel += 1
                else:
                    aantal_niet += 1
        self.__x = ["aantal hypothetical protein", "anders"]
        self.__y = [aantal_wel, aantal_niet]

    def setgegevens_strands(self):
        """
        Berekent de gegevens voor de + en - strands
        """
        self.__plus_grote = 0
        self.__min_grote = 0
        try:
            for lijst in self.__gff_totale_inhoud:
                if lijst[6] == "+":
                    self.__plus_grote += 1
                elif lijst[6] == "-":
                    self.__min_grote += 1
        except IndexError:  # in case er geen lijst[6] is
            print("list out of index!")
        self.__grote = [self.__plus_grote, self.__min_grote]
        totaal = self.__grote[0] + self.__grote[1]
        # zodat de cirkel vol is
        self.__grote[0] = (self.__grote[0] / totaal)
        self.__grote[1] = (self.__grote[1] / totaal)

    def setLengte(self):
        """
         Kijkt hoevaak elk lengte voorkomt en voegt het toe aan een dic.
        """
        self.__verschillende_lengtes = {}
        for gene in self.__gff_gene:
           lengte = int(gene[4]) - int(gene[3]) # berekent de lengte
           # voegt alleen unieke lengtes toe of voeg 1 toe aan de waarde
           if lengte not in self.__verschillende_lengtes:
               self.__verschillende_lengtes[lengte] = 1
           else:
               self.__verschillende_lengtes[lengte] += 1

    def maakGrafiek_lengte(self):
        """
        Laat zien hoevaak elke lengte voor komt van de genen
        """
        scherm = Tk.Toplevel() # zorgt ervoor dat de tabel op een
        # nieuw GUI scherm komt
        scherm.title("verschillende lengtes van genen")
        scrollbar = Tk.Scrollbar(scherm) # maakt een scroll bar
        scrollbar.pack(side=Tk.RIGHT, fill = Tk.Y)
        gegevens = Tk.Listbox(scherm, yscrollcommand = scrollbar.set
                              ,width=20)
        # voegt de gegevens toe aan de scroll bar
        for k,v in self.__verschillende_lengtes.items():
            zin = "lengte:", k ,"=", v, "x"
            gegevens.insert(Tk.END, zin)
        gegevens.pack(side=Tk.LEFT,fill=Tk.BOTH)
        scrollbar.config(command=gegevens.yview)
        # return knop
        knop = Tk.Button(text="return", master=scherm,
                         command=scherm.destroy)
        knop.pack()

    def tabel_patronen(self):
        """
        maakt de tabel aan met de gevonden patronen
        """
        tabelscherm = Tk.Toplevel()  # zorgt ervoor dat de tabel op een
        # nieuw GUI scherm komt
        tabelscherm.title("gevonden patronen")
        aantal_rijen = len(self.__tabel)
        aantal_kolommen = len(self.__tabel[0])
        # loop wat zorgt voor een juist aantal rijen en kolommen
        # en voert de inhoud toe
        for rij in range(aantal_rijen):
            for kolom in range(aantal_kolommen):
                tabel = Tk.Entry(tabelscherm, fg="black", width=50)
                tabel.grid(row=rij, column=kolom)
                tabel.insert(Tk.END, self.__tabel[rij][kolom])
        knop = Tk.Button(text="return", master=tabelscherm,
                         command=tabelscherm.destroy)
        knop.grid()

    def grafiek_type(self):
        """
        Maakt de grafiek van verschillende types in de GUI
        """
        grafiek_type = Tk.Toplevel()  # Zorgt ervoor dat de grafiek
        # op een ander scherm komt
        keys = []  # x-as waarde
        data = []  # y-as waarde
        for k, v in self.__labels.items():
            data.append(v)
            keys.append(k)

        f = Figure(figsize=(10, 8), dpi=100)
        ax = f.add_subplot(111)
        width = .8
        # maakt de grafiek
        rects1 = ax.bar(keys, data, width)
        figuur = FigureCanvasTkAgg(f, master=grafiek_type)
        figuur.draw()
        figuur.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        quit_button = Tk.Button(heigh=2, width=10, text="return",
                                command=grafiek_type.destroy,
                                master=grafiek_type)
        quit_button.pack()

    def grafiek_hypothetical_proteins(self):
        """
        Maakt de grafiek van hypothetical proteins
        """
        grafiek_hypothetical_proteins = Tk.Toplevel()  # zorgt ervoor
        # dat de grafiek op een nieuw venster komt te staan
        # maakt de grafiek
        f = Figure(figsize=(10, 8), dpi=100)
        ax = f.add_subplot(111)
        width = .8
        rects1 = ax.bar(self.__x, self.__y, width)
        figuur = FigureCanvasTkAgg(f, master=grafiek_hypothetical_proteins)
        figuur.draw()
        figuur.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        # zorgt voor een return button
        quit_button = Tk.Button(heigh=2, width=10, text="return",
                                command=grafiek_hypothetical_proteins.destroy,
                                master=grafiek_hypothetical_proteins)
        quit_button.pack()

    def maak_grafiek_strands(self):
        """
        Maakt de grafiek met de hoeveelheid +/- strands
        """
        grafiek_strands = Tk.Toplevel()  # zorgt ervoor de de grafiek op
        # een nieuw venster komt te staan
        upper_frame = Tk.Frame(grafiek_strands)
        bottom_frame = Tk.Frame(grafiek_strands)
        upper_frame.pack()
        bottom_frame.pack()
        label = Tk.Label(upper_frame, text='percentage +/-strands')
        # legenda met betekenis van de kleuren
        legenda = Tk.Label(bottom_frame, text=("legenda:\n "
                                               "- blauw: +", self.__plus_grote,
                                               "\n-paars: -",
                                               self.__min_grote))
        label.pack()
        legenda.pack()
        pie = Tk.Canvas(upper_frame, width=154, height=154)
        pie.pack()
        # maakt de + deel
        pie.create_arc((2, 2, 152, 152), fill="#00BFFF", start=0,
                       extent=(360 * self.__grote[0]))
        # maakt de - deel
        pie.create_arc((2, 2, 152, 152), fill="#8A2BE2",
                       start=(360 * self.__grote[0]),
                       extent=(360 * self.__grote[1]))
        # maakt return button
        quit_button = Tk.Button(heigh=2, width=10, text="return",
                                command=grafiek_strands.destroy,
                                master=bottom_frame)
        quit_button.pack()

    def maak_grafiek_percentage(self):
        """
        Maakt de grafiek van GC en AT %
        """
        grafiek_percentage = Tk.Toplevel()  # zorgt ervoor de de grafiek op
        # een nieuw venster komt te staan
        upper_frame = Tk.Frame(grafiek_percentage)  # frame voor grafiek
        bottom_frame = Tk.Frame(grafiek_percentage)  # frame voor legenda
        upper_frame.pack()
        bottom_frame.pack()
        label = Tk.Label(upper_frame, text='percentage AT/GC %')
        # legenda met betekenis en de percentage
        legenda = Tk.Label(bottom_frame, text=(
            "legenda:\n - blauw: GC:", self.__percentage[0],
            "\n- paars: AT: - ", self.__percentage[1]))
        label.pack()
        legenda.pack()
        # maakt de circel aan
        pie = Tk.Canvas(upper_frame, width=154, height=154)
        pie.pack()
        pie.create_arc((2, 2, 152, 152), fill="#00BFFF", start=0,
                       extent=(360 * self.__gc) / 100)
        pie.create_arc((2, 2, 152, 152), fill="#8A2BE2",
                       start=(360 * self.__gc) / 100,
                       extent=(360 * self.__at) / 100)
        # return button
        quit_button = Tk.Button(heigh=2, width=10, text="return",
                                command=grafiek_percentage.destroy,
                                master=bottom_frame)
        quit_button.pack()


class Patronen():
    def __init__(self, cds, gffb_gene, gffb_cds):
        self.__gff_cds = cds
        self.__gffb_cds = gffb_cds
        self.__gffb_gene = gffb_gene
        self.vindpatronen()
        self.vindprotein_id()
        self.vindproduct()
        self.vindstrand()

    def vindpatronen(self):
        """
        zoekt naar de patronen en zet ze in een 2d lijst
        """
        self.__gevonden_patronen = []
        aantal = 1
        for key, value in self.__gffb_cds.items():
            patroon = []
            for i in value:  # gaat alle values langs todat die /translation=
                # heeft gevonden
                if i.startswith("/translation"):
                    x = bool(re.match(  # patroon 1
                        ".*[ST]G[LIVMFYW]{3}[GN][A-Z]{2}T[LIVM][A-Z]T[A-Z]{2}H.*",
                        i))
                    y = bool(  # patroon 2
                        re.match(".*T[A-Z]{2}[GC][NQ]SGS[A-Z][LIVM][FY].*", i))
                    if x == True or y == True:
                        patroon.append(aantal)
                        i = i.replace("/translation=", "")
                        patroon.append(i)
            # voegt de lijst toe aan de grote lijst
            if patroon != []:
                self.__gevonden_patronen.append(patroon)
                aantal += 1

    def vindprotein_id(self):
        """
        Zoekt bij alle gevonden patronen de correcte protein_id
        """
        for patroon in self.__gevonden_patronen:
            value = self.__gffb_cds[patroon[0]]  # kijkt naar de id van de
            # gevonden patroon en zoekt op dat id in de dic. CDS
            for i in value:
                if i.startswith("/protein_id="):
                    i = i.replace("/protein_id=", "")  # zorgt dat je alleen
                    # de protein_id leest
                    patroon.append(i)

    def vindproduct(self):
        """
        zoekt naar het product voor de gevonden patronen
        """
        for patroon in self.__gevonden_patronen:
            value = self.__gffb_cds[patroon[0]]  # kijkt naar de id van de
            # gevonden patroon en zoekt op dat id in de dic. CDS
            for i in value:
                if i.startswith("/product="):
                    i = i.replace("/product=", "")
                    patroon.append(i)

    def vindstrand(self):
        """
        Zoekt naar de strand die hoort bij de gevonden patronen
        """
        for patroon in self.__gevonden_patronen:
            key = patroon[0]
            value = self.__gff_cds[key]
            strand = value[6]
            patroon.append(strand)

    def getpatronen(self):
        """

        :return: self.__gevonden_patronen - 2d lijst met alle gevonden patronen
        [[1,"MKGKFL..", "YP_498609.1", "chromosomal ...","+"], [2, ..] n]
        """
        return self.__gevonden_patronen


def sequentie_uit_bestand3(bestandnaam):
    """
    zorgt dat de sequentie achter elkaar in een string komt te staan
    ( oorspronkelijk hoort deze code in de class inhoud_gbff, maar doordat
    het script veel langzamer werd als het daarin stond heb ik het maar
    buiten de class neergezet )
    :param bestandnaam: bestandnaam - naam van bestand 3
    :return: sequentie - string met de sequentie
    """
    bestand = open(bestandnaam.name, "r")
    sequentie = ""
    for regel in bestand:
        for index in regel:
            if index == "a" or index == "t" or index == "c" \
                    or index == "g":  # als de index een nucleotide is wordt
                # het toegevoegt aan de sequentie
                sequentie += index
    bestand.close()
    return sequentie


def main():
    bestandsnaam_gff = "gff_bestand.txt"
    bestandsnaam_gbff = "gbff.txt"

    gff = Gff_bestand(bestandsnaam_gff)
    gff_cds, gff_gene, gff_trna, ggf_exon, gff_anders, \
    gff_totale_inhoud = gff.getlijsten()
    gffb = Gffb_bestand(bestandsnaam_gbff)
    gffb_gene, gffb_cds, gffb_rrna = gffb.getdiconaries()
    sequentie_bestand = gffb.getsequentie()
    sequentie = sequentie_uit_bestand3(sequentie_bestand)

    p = Patronen(gff_cds, gffb_gene, gffb_cds)
    patronen = p.getpatronen()

    g = Gui_maken(patronen, gff_totale_inhoud, gffb_cds, sequentie,
                  gff_gene)


main()


# to do:
# - meer RE toepassen
# - als het lukt overerving
# - pep8 regels controleren
