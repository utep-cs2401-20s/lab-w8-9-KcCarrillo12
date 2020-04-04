import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

public class AminoAcidLLTester {

    @Test
    public void codonCompare(){
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence("AAGAUGCCU");
        AminoAcidLL.sort(sol);
        int expected = 0;
        assertEquals(expected, sol.codonCompare(sol));
    }

    @Test
    public void codonCompare2() {
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence("CAUCAC");
        AminoAcidLL.sort(sol);
        int expected = 0;
        assertEquals(expected, sol.codonCompare(sol));
    }

    @Test
    public void aminoAcidList() {
        String RNASequence = "AAGAUGCCU";
        AminoAcidLL list = AminoAcidLL.createFromRNASequence(RNASequence);
        char[] testArr = {'K','M','P'};
        char[] result = list.aminoAcidList();
        assertArrayEquals(testArr, result);
    }

    @Test
    public void aminoAcidList2() {
        String RNASequence = "CAUCAC";
        AminoAcidLL list = AminoAcidLL.createFromRNASequence(RNASequence);
        char[] testArr = {'H'};
        char[] result = list.aminoAcidList();
        assertArrayEquals(testArr, result);
    }

    @Test
    public void isSorted(){
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence("AAGAUGCCGU");
        assertEquals(true, sol.isSorted());
    }

    @Test
    public void isSorted2(){
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence("CAUCAC");
        assertEquals(true, sol.isSorted());
    }

    @Test
    public void createFromRNASequence() {
        String RNASequence = "AAGAUGCCGU";
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence(RNASequence);
    }

    @Test
    public void createFromRNASequence2() {
        String RNASequence = "CAUCAC";
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence(RNASequence);
    }

    @Test
    public void sort(){
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence("AAGAUGCCGU");
        sol = AminoAcidLL.sort(sol);
        assertEquals(true, sol.isSorted());
    }

    @Test
    public void sort2(){
        AminoAcidLL sol = AminoAcidLL.createFromRNASequence("CAUCAC");
        sol = AminoAcidLL.sort(sol);
        assertEquals(true, sol.isSorted());
    }
}