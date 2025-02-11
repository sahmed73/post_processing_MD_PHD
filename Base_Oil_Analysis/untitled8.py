class Solution(object):
    def addTwoNumbers(self, l1, l2):
        """
        :type l1: ListNode
        :type l2: ListNode
        :rtype: ListNode
        """
        l3 = []
        hand = 0
        for d1,d2 in zip(l1[::-1],l2[::-1]):
            dsum = d1+d2+hand
            l3.append(dsum%10)
            hand = dsum//10
        l3.append(hand)
        return l3[::-1]

sol = Solution()
l1 = [9,9,9,9,9,9,9]
l2 = [0,0,0,9,9,9,9]
sol.addTwoNumbers(l1, l2)