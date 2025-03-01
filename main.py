from flask import Flask, request, render_template_string

app = Flask(__name__)

def calcule(x):
    if isinstance(x, int) and x >= 0:  # Vérifie que x est un entier positif
        F = 1
        for i in range(1, x + 1):
            F *= i
        return F
    return None  # Retourne None si x n'est pas valide

@app.route("/", methods=["GET", "POST"])
def index():
    result = ""
    if request.method == "POST":
        try:
            x = int(request.form["x"])
            y = calcule(x)
            if y is not None:
                result = f"<p>Le factoriel de {x} est : {y}</p>"
            else:
                result = "<p>Veuillez entrer un entier positif.</p>"
        except ValueError:
            result = "<p>Entrée invalide. Veuillez entrer un nombre entier.</p>"

    return render_template_string("""
    <html>
    <head>
        <title>Calcul de Factorielle</title>
    </head>
    <body>
        <h2>Calcul de Factorielle</h2>
        <form method="POST">
            <label for="x">x :</label>
            <input type="text" name="x" required><br><br>
            <button type="submit">Calculer</button>
        </form>
        <br>
        {{ result | safe }} 
    </body>
    </html>""", result=result)

if __name__ == "__main__":
    app.run(debug=True)