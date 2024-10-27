# Non-Obtuse-Triangulation-PSLG

Using CGAL library for triangulation proccess of a PSLG targeting to non-obtuse angles

## Ομάδα 38

- ### ΦΟΥΤΡΗΣ ΒΑΣΙΛΕΙΟΣ 1115202000231

- ### ΠΑΠΑΔΗΜΗΤΡΙΟΥ ΜΑΡΙΝΑ 1115202100136

## Λειτουργίες

### `void parseInput(const nlohmann::json& input, std::vector<Point>& points, std::vector<std::pair<int, int>>& constraints)`

Αυτή η συνάρτηση αναλύει τα δεδομένα εισόδου από ένα JSON αντικείμενο. Ελέγχει την ύπαρξη των κλειδιών `points_x` και `points_y`, και εάν η ποσότητα των σημείων που δίνονται συμφωνεί με τον αριθμό που έχει καθοριστεί. Επίσης, διαβάζει τυχόν επιπλέον περιορισμούς από το JSON.

### `double approximate_angle_2D(const Point& p1, const Point& p2, const Point& p3)`

Υπολογίζει τη γωνία μεταξύ τριών σημείων σε 2D χρησιμοποιώντας το νόμο των συνημίτονων. Επιστρέφει τη γωνία σε μοίρες.

### `bool isObtuseAngle(const Point& a, const Point& b, const Point& c)`

Ελέγχει αν η γωνία που σχηματίζεται από τρία σημεία είναι οξεία. Χρησιμοποιεί τον τριγωνικό ανισότητα για να προσδιορίσει την οξύτητα.

### `bool isObtuse(CDT::Face_handle face)`

Ελέγχει αν κάποιο από τα τρία τρίγωνα μιας δεδομένης όψης είναι οξεία. Επιστρέφει `true` αν υπάρχει οξεία γωνία.

### `Point obtuseAngleProjection(const Point& obtuse_vertex, const Point& p1, const Point& p2)`

Προβάλει την οξεία κορυφή σε μια γραμμή μεταξύ των δύο άλλων κορυφών του τριγώνου. Χρησιμοποιείται όταν ανιχνεύεται οξεία γωνία.

### `Point insertSteinerAtMedian(const Point& p1, const Point& p2)`

Εισάγει ένα Steiner σημείο στο μέσο της γραμμής μεταξύ δύο σημείων, δηλαδή στην μέση του μεγαλύτερου άκρου.

### `Point insertSteinerAtCircumcenter(const Point& p1, const Point& p2, const Point& p3)`

Εισάγει ένα Steiner σημείο στο κυκλικό κέντρο του τριγώνου που σχηματίζεται από τρία σημεία.

### `Point insertSteinerForPolygon(const std::vector<Point>& polygon_points)`

Εισάγει ένα Steiner σημείο στο κέντρο βάρους ενός πολυγώνου, υπολογίζοντας το μέσο όρο των συντεταγμένων των σημείων του πολυγώνου.

### `void insertSteinerPoint(CDT& cdt, CDT::Face_handle face)`

Εισάγει ένα Steiner σημείο σε μια δεδομένη όψη ελέγχοντας πρώτα αν η όψη περιέχει οξεία γωνία και επιλέγοντας τη στρατηγική εισαγωγής ανάλογα.

### `void applyPolygonalTreatment(CDT& cdt, CDT::Face_handle face1, CDT::Face_handle face2)`

Εφαρμόζει μια λυση σε πολυγωνικές περιοχές που σχηματίζονται από δύο γειτονικές όψεις, εισάγοντας ένα Steiner σημείο στο κέντρο βάρους των κορυφών τους.

### `Face_handle getAdjacentFace(Face_handle face, int edge_index)`

Επιστρέφει την γειτονική όψη που βρίσκεται απέναντι από μια συγκεκριμένη ακμή στην δεδομένη όψη.

### `CDT::Face_handle findObtuseAdjacentFace(CDT& cdt, CDT::Face_handle face)`

Αναζητά μια γειτονική όψη που είναι οξεία, επιστρέφοντας τη γειτονική όψη εάν βρεθεί.

### `bool isEdgeFlipBeneficial(CDT& cdt, CDT::Face_handle& face, int edge_index)`

Ελέγχει αν η ανατροπή μιας ακμής θα είναι επωφελής, υπολογίζοντας τις γωνίες των νέων τριγώνων που θα σχηματιστούν.

### `bool applyEdgeFlip(CDT& cdt, CDT::Face_handle& face)`

Εφαρμόζει την ανατροπή των ακμών στην όψη, εάν κρίνεται επωφελής, και επιστρέφει αν η ανατροπή ήταν επιτυχής.

### `void writeOutput(const CDT& cdt, const std::vector<Point>& steinerPoints, nlohmann::json& output)`

Γράφει τα αποτελέσματα της τριγωνοποίησης σε μορφή JSON, συμπεριλαμβανομένων των όψεων και των Steiner σημείων.

## Σημείωση

Αυτή η εφαρμογή απαιτεί τη βιβλιοθήκη CGAL και nlohmann/json για την ανάλυση και την έξοδο των δεδομένων. Βεβαιωθείτε ότι έχετε εγκαταστήσει αυτές τις βιβλιοθήκες πριν από την εκτέλεση του προγράμματος.

**To run manually:**
-del dir built/mkdir built
-cd built
-cmake ..
-make
-./triangulation_program

**For the run.sh :**
chmod +x run.sh

**For the cleanup.sh :**
chmod +x cleanup.sh

Before commitng changes run the cleanup.sh

[Non-Obtuse Triangulation PSLG](https://github.com/HackeRinaa/Non-Obtuse-Triangulation-PSLG)
