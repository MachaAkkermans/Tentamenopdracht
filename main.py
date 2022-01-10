import matplotlib.pyplot as plt
import re

class Grafieken():

    def __init__(self):
        self.__grafieksoort = input("Welk grafiek wil je zien?\n\n"
                                    "Kies uit:\n- aantal sequenties"
                                    "+/- strand\n - ")


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
        self.__title = "GC en AT %"
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

    def __init__(self, inhoud_gff):
        self.__inhoud_gff = inhoud_gff
        self.setgegevens()

    def setgegevens(self):
        self.__labels = {}
        for lijst in self.__inhoud_gff:
            if lijst[2] not in self.__labels:
                self.__labels[lijst[2]] = 1
            else:
                self.__labels[lijst[2]] += 1
        print(self.__labels)

    def maakgrafiek(self):
        plt.bar(range(len(self.__labels)), self.__labels.values())
        plt.title("grafiek met hoevaak elke type voorkomt")
        plt.xlabel("Type")
        plt.ylabel("aantal")
        plt.grid()
        plt.show()



def grafieken(inhoud_gff,sequentie_genoom):
    grafieksoort = input("Welk grafiek wil je zien?\n\nKies uit:\n"
                         "- aantal sequenties +/- strand\n- types\n"
                         "- GC%\n")
    if grafieksoort == "aantal sequenties +/- strand":
        g = Taartdiagram(inhoud_gff)
        g.setgegevens_strands()
        print (g.maakgrafiek())
    elif grafieksoort == "types":
        g = Staafdiagram(inhoud_gff)
        print (g.maakgrafiek())
    elif grafieksoort == "GC%":
        print("jodiejo")
        g = Taartdiagram(sequentie_genoom)
        g.setgegevens_gc()
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

#def informatielezen_gbff(bestandsnaam):

def gbff_bestand_opdelen(bestandsnaam_gbff):
    #bestand1 = open("bestand1","w")
    #bestand2 = open("gene+cds.txt","w")
    #bestand3 = open("sequentie.txt","w")
#
    #bestand = open(bestandsnaam_gbff,"r")
    #regel_nummer = 0
    #for regel in bestand:
    #    if regel_nummer >= 0 and regel_nummer <= 45:
    #        bestand1.write(regel)
    #    elif regel_nummer <= 49099 and regel_nummer >= 46:
    #        bestand2.write(regel)
    #    elif regel_nummer >= 49100:
    #        bestand3.write(regel)
    #    regel_nummer += 1
    #bestand1.close()
    #bestand2.close()
    #bestand3.close()
    #return bestand2, bestand3
#
def #sequentie_uit_bestand(sequentie_bestand):
    #bestand = open("sequentie.txt", "r")
    #sequentie = ""
    #for regel in bestand:
    #    for index in regel:
    #        if index == "a" or index == "t" or index == "c" or index == "g":
    #            sequentie += index
    #return sequentie
#

def main():
    #bestandsnaam = "test.txt"
    bestandsnaam_gff = "gff_bestand.txt"
    bestandsnaam_gbff = "gbff.txt"
    inhoud_gff = informatieinlezen_gff(bestandsnaam_gff)
    gene_cds_bestand, sequentie_bestand = gbff_bestand_opdelen(bestandsnaam_gbff)
    sequentie_genoom = sequentie_uit_bestand(sequentie_bestand)
    grafieken(inhoud_gff,sequentie_genoom)


main()
