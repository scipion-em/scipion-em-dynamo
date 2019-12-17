from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from tomo.protocols import ProtImportSubTomograms
from tomo.tests import DataSet

from dynamo.protocols import DynamoSubTomoMRA


class TestSubTomogramsAlignment(BaseTest):
    """ This class check if the protocol to import sub tomograms works
    properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runAlignment(self, randomInitialization=True, numberOfReferences=2,
                      numberOfIters=3, angularSampling=15):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)

        alignment = self.newProtocol(DynamoSubTomoMRA,
                                      inputVolumes=protImport.outputSubTomograms)
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.outputSubtomograms,
                             "There was a problem with SetOfSubtomograms output")
        self.assertIsNotNone(alignment.outputClassesSubtomo,
                             "There was a problem with outputClassesSubtomo output")
        return alignment

    def test_outputMltomo(self):
        dynamoAlignment = self._runAlignment()
        outputSubtomos = getattr(dynamoAlignment, 'outputSubtomograms')
        outputClasses = getattr(dynamoAlignment, 'outputClassesSubtomo')
        self.assertTrue(outputSubtomos)
        self.assertTrue(outputSubtomos.getFirstItem().hasTransform())
        self.assertTrue(outputClasses)
        self.assertTrue(outputClasses.hasRepresentatives())
        return dynamoAlignment
