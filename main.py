from tkinter import messagebox
import matplotlib.pyplot as plt
import re
import tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class Grafieken():

    def __init__(self):
        self.__grafieksoort = input("Welk grafiek wil je zien?\n\n"
                                    "Kies uit:\n- aantal sequenties"
                                    "+/- strand\n - ")


    def labels(x,y):
        for i in range(len(x)):
            plt.text(i, y[i], y[i])


class Taartdiagram(Grafieken):

    def __init__(self, gegevens):
        self.__gegevens = gegevens
        #self.setgegevens_strands()


    def setgegevens_strands(self):
        self.__title = "Aantal sequenties op +/- strand"
        self.__labels = "+","-"
        self.__grote = [0,0]
        for lijst in self.__gegevens:
            if lijst[6] == "+":
                self.__grote[0] += 1
            elif lijst[6] == "-":
                self.__grote[1] += 1

    def setgegevens_gc(self):
        self.__title = "GC/AT %"
        self.__labels = "GC%", "AT%"
        gc = round((self.__gegevens.count("c") + self.__gegevens.count("g")) / len(self.__gegevens),2) *100
        at = round((self.__gegevens.count("a") + self.__gegevens.count("t")) / len(self.__gegevens),2) *100
        self.__grote = [gc,at]

    def maakgrafiek(self):
        figuur, vormgeving = plt.subplots()
        vormgeving.pie(self.__grote, labels=self.__labels)
        plt.title(self.__title)
        plt.text(0,0.6, self.__grote[0])
        plt.text(0,-0.6, self.__grote[1])
        plt.show()


class Staafdiagram(Grafieken):

    def __init__(self, gegevens):
        self.__gegevens = gegevens

    def labels(self):
        for i in range(len(self.__x)):
            plt.text(i, self.__y[i], self.__y[i])

    def setgegevenstypes(self):
        self.__labels = {}
        for lijst in self.__gegevens:
            if lijst[2] not in self.__labels:
                self.__labels[lijst[2]] = 1
            else:
                self.__labels[lijst[2]] += 1
        self.__x = []
        self.__y = []
        for key in self.__labels.keys():
            self.__x.append(key)
        for value in self.__labels.values():
            self.__y.append(value)
        self.labels()


    def setgegevenshypothetical_proteins(self):
        aantal = 0
        for k,v in self.__gegevens.items():
            for item in v:
                if item == "/note=\"conserved hypothetical protein\"":
                    aantal += 1
        self.__x = ["aantal hypothetical protein"]
        self.__y = [aantal]
        self.labels()

    def maakgrafiek(self):
        plt.bar(self.__x,self.__y)
        plt.title("grafiek met hoevaak elke type voorkomt")
        plt.xlabel("Type")
        plt.ylabel("aantal")
        plt.grid()
        plt.show()


class gui_maken():
    def __init__(self,patronen,inhoud_gff):
        self.__patronen = patronen
        self.__inhoud_gff = inhoud_gff
        self.setgegevenstypes()
        self.__tabel = []
        self.settabel()
        self.gui_startscherm()
        #self.grafiek_type()

    def setgegevenstypes(self):
        self.__labels = {}
        for lijst in self.__inhoud_gff:
            if lijst[2] not in self.__labels:
                self.__labels[lijst[2]] = 1
            else:
                self.__labels[lijst[2]] += 1
        self.__x = []
        self.__y = []
        for key in self.__labels.keys():
            self.__x.append(key)
        for value in self.__labels.values():
            self.__y.append(value)

        self.labelsgrafiek()

    def labelsgrafiek(self):
        for i in range(len(self.__x)):
            plt.text(i, self.__y[i], self.__y[i])

    def settabel(self):
        self.__tabel.append(["nummer","eiwit","protein_id","product","+/- strand"])
        for patroon in self.__patronen:
            self.__tabel.append(patroon)


    def gui_startscherm(self):
        startscherm = Tk.Tk()
        self.__knop = Tk.Button(heigh=2, master=startscherm, text="klik hier voor tabel patronen",command=self.tabel_patronen)
        self.__knop2 = Tk.Button(heigh=2, master=startscherm, text="klik hier voor grafiek met verschillende type ",command=self.grafiek_type)
        self.__knop3 = Tk.Button(heigh=2, master=startscherm,text="klik hier voor grafiek met verschillende type ",command=self.grafiek_type)
        self.__knop.pack()
        self.__knop2.pack()
        self.__knop3.pack()
        Tk.mainloop()

    def tabel_patronen(self):
        tabelscherm =Tk.Toplevel()
        tabelscherm.title("gevonden patronen")
        aantal_rijen = len(self.__tabel)
        aantal_kolommen = len(self.__tabel[0])
        for rij in range(aantal_rijen):
            for kolom in range(aantal_kolommen):
                tabel =Tk.Entry(tabelscherm,fg="black",width=50)
                tabel.grid(row=rij,column=kolom)
                tabel.insert(Tk.END, self.__tabel[rij][kolom])
        knop = Tk.Button(text="return", master=tabelscherm, command=tabelscherm.destroy)
        knop.grid()

    def grafiek_type(self):
        grafiek_type =Tk.Toplevel()

        f = Figure(figsize=(10, 8), dpi=100)
        ax = f.add_subplot(111)

        keys = []
        data = []
        for k, v in self.__labels.items():
            data.append(v)
            keys.append(k)

        width = .8

        rects1 = ax.bar(keys, data, width)

        figuur = FigureCanvasTkAgg(f, master=grafiek_type)
        figuur.draw()
        figuur.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        quit_button = Tk.Button(heigh=2, width=10, text="return",
                                command=grafiek_type.destroy, master=grafiek_type)
        quit_button.pack()


class Patronen():
    def __init__(self,cds,inhoud_gff):
        self.__cds = cds
        self.__inhoud_gff = inhoud_gff
        del self.__inhoud_gff[0] #de region regio heb je niet nodig
        self.vindpatronen()
        self.vindprotein_id()
        self.vindproduct()
        self.vindstrand()


    def vindpatronen(self):
        self.__gevonden_patronen = []
        aantal = 1
        for key, value in self.__cds.items():
            patroon = []
            for i in value:
                if i.startswith("/translation"):
                    x = bool(re.match(
                        ".*[ST]G[LIVMFYW]{3}[GN][A-Z]{2}T[LIVM][A-Z]T[A-Z]{2}H.*",
                        i))
                    y = bool(
                        re.match(".*T[A-Z]{2}[GC][NQ]SGS[A-Z][LIVM][FY].*", i))
                    if x == True or y == True:
                        patroon.append(aantal)
                        i = i.replace("/", "")
                        patroon.append(i)
            if patroon != []:
                self.__gevonden_patronen.append(patroon)
                aantal += 1

    def vindprotein_id(self):
        for patroon in self.__gevonden_patronen:
            value = self.__cds[patroon[0]]
            for i in value:
                if i.startswith("/protein_id="):
                    i = i.replace("/", "")
                    patroon.append(i)

    def vindproduct(self):
        for patroon in self.__gevonden_patronen:
            value = self.__cds[patroon[0]]
            for i in value:
                if i.startswith("/product="):
                    i = i.replace("/","")
                    patroon.append(i)

    def vindstrand(self):
        for patroon in self.__gevonden_patronen:
            key = patroon[0]
            key = key * 2 #in het gff bestand staat het als - gene
                                                          # - cds
                                                          # - gene
                          # dus keer 2 zodat je de goed cds hebt
            value = self.__inhoud_gff[key]
            strand = value[6]
            patroon.append(strand)

    def getpatronen(self):
        return self.__gevonden_patronen


def grafieken(inhoud_gff,sequentie_genoom,cds):
    grafieksoort = input("Welk grafiek wil je zien?\n\nKies uit:\n"
                         "- aantal sequenties +/- strand\n- types\n"
                         "- GC%\n- aantal hypothetical proteins\n")
    if grafieksoort == "aantal sequenties +/- strand":
        g = Taartdiagram(inhoud_gff)
        g.setgegevens_strands()
        print (g.maakgrafiek())
    elif grafieksoort == "types":
        g = Staafdiagram(inhoud_gff)
        g.setgegevenstypes()
        print (g.maakgrafiek())
    elif grafieksoort == "GC%":
        g = Taartdiagram(sequentie_genoom)
        g.setgegevens_gc()
        print (g.maakgrafiek())
    elif grafieksoort == "aantal hypothetical proteins":
        g = Staafdiagram(cds)
        g.setgegevenshypothetical_proteins()
        print (g.maakgrafiek())

    else:
        print("Verzoek niet gevonden! probeer het opnieuw")
        grafieken(inhoud_gff)


def informatieinlezen_gff(bestandsnaam):
    bestand = open(bestandsnaam)
    inhoud_gff = []
    for regel in bestand:
        regel = regel.replace("\n","") # Haalt de enter eruit
        lijst = regel.split("\t")      # splits de lijst op elke tab
        inhoud_gff.append(lijst)
    bestand.close()
    return inhoud_gff


def gbff_bestand_opdelen(bestandsnaam_gbff):
    bestand1 = open("bestand1","w")
    bestand2 = open("gene+cds.txt","w")
    bestand3 = open("sequentie.txt","w")

    bestand = open(bestandsnaam_gbff,"r")
    regel_nummer = 0
    for regel in bestand:
        if regel_nummer >= 0 and regel_nummer <= 45:
            bestand1.write(regel)
        elif regel_nummer <= 49099 and regel_nummer >= 46:
            bestand2.write(regel)
        elif regel_nummer >= 49100:
            bestand3.write(regel)
        regel_nummer += 1
    bestand1.close()
    bestand2.close()
    bestand3.close()
    return bestand2, bestand3


def sequentie_uit_bestand(sequentie_bestand):
    bestand = open("sequentie.txt", "r")
    sequentie = ""
    for regel in bestand:
        for index in regel:
            if index == "a" or index == "t" or index == "c" or index == "g":
                sequentie += index
    return sequentie


def informatie_uit_gene_cds_bestand(gene_cds_bestand):
    """

    :return: lijst_met_gene_cds - 2d lijst met inhoudt van bestand gene+cds.txt
    - [["gene 517..1878", /gene="dnaA", ...] [CDS 517..1878, /gene="dnaA",...],
      ["gene ..., n]]
    """
    bestand = open(gene_cds_bestand.name)
    #bestand = open("test.txt")
    lijst_met_gene_cds = []
    per_cds_of_gene = []
    per_kopje_in_bestand = ""
    aantal = 0
    for regel in bestand:
        while regel.startswith(" "):  #Haalt alle spaties weg dat bij het
            regel = regel[1:]         # begin van de regel staat
        regel = regel.replace("\n","")
        # zorgt ervoor dat elke keer wanneer "gene" of "CDS"staat er een nieuw
        # lijst wordt aangemaakt en de oude wordt toegevoegt aan de 2d lijst
        if regel.startswith("gene") or regel.startswith("CDS"):
            per_cds_of_gene.append(per_kopje_in_bestand)
            lijst_met_gene_cds.append(per_cds_of_gene)
            per_cds_of_gene = []
            per_cds_of_gene.append(regel)
            per_kopje_in_bestand = ""
        else:
            # Zorgt ervoor dat kopjes die bestaan uit meerdere regels
            # samen als een regel worden toegevoegt in een lijst
            if regel.startswith("/"):
                if per_kopje_in_bestand != "":
                    per_cds_of_gene.append(per_kopje_in_bestand)
                    per_kopje_in_bestand = ""
                per_kopje_in_bestand += regel
            else:
                per_kopje_in_bestand += regel
    per_cds_of_gene.append(per_kopje_in_bestand)
    lijst_met_gene_cds.append(per_cds_of_gene)
   #    aantal += 1
   #     if aantal == 50:
   #         break

    del lijst_met_gene_cds[0]
    #print (lijst_met_gene_cds)
    return lijst_met_gene_cds


def maken_dictonary(lijst_met_gene_cds):
    aantal = 0
    gene = {}
    cds = {}
    id_nummer = 1
    for lijst in lijst_met_gene_cds:
        if aantal%2 == 0:
            gene[id_nummer] = lijst
            #id_nummer += 1
        else:
            cds[id_nummer] = lijst
            id_nummer += 1
        aantal += 1
    return gene, cds


def main():
    bestandsnaam_gff = "gff_bestand.txt"
    bestandsnaam_gbff = "gbff.txt"

    inhoud_gff = informatieinlezen_gff(bestandsnaam_gff)
    gene_cds_bestand, sequentie_bestand = gbff_bestand_opdelen(bestandsnaam_gbff)

    sequentie_genoom = sequentie_uit_bestand(sequentie_bestand)
    lijst_met_gene_cds = informatie_uit_gene_cds_bestand(gene_cds_bestand)
    gene, cds = maken_dictonary(lijst_met_gene_cds)

    #grafieken(inhoud_gff,sequentie_genoom,cds)
    #gevonden_patronen = patronen(cds)
    p = Patronen(cds,inhoud_gff)
    patronen = p.getpatronen()

    g = gui_maken(patronen,inhoud_gff)

main()
