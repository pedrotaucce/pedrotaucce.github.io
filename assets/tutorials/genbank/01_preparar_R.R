# =====================================================================
# Prática 3 - Aquisição de sequências no GenBank
# Script 1: preparar o R (executar uma única vez em cada computador)
# =====================================================================

# Instala os dois pacotes usados na prática:
#   ape     - leitura, escrita e análise de sequências (Paradis & Schliep 2019)
#   rentrez - acesso às bases do NCBI pelo R (Winter 2017)
install.packages(c("ape", "rentrez"))

# Confere a instalação: library() carrega cada pacote. Se não aparecer
# nenhuma mensagem de erro (como "there is no package called 'ape'"),
# está tudo certo.
library(ape)
library(rentrez)

# Versões instaladas: devem aparecer como [1] '5.8.1', por exemplo.
# Anote-as: são parte da documentação do seu trabalho.
packageVersion("ape")
packageVersion("rentrez")

# Pasta de trabalho atual: é aqui que os arquivos serão gravados
getwd()
